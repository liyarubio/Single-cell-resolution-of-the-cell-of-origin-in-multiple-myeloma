library(ArchR)
library(tidyr)
library(parallel)
library(ChIPseeker)
library(Cairo)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggpubr)
## ggplot2 3.4.0

setwd("/media/liyaru/LYR/MM2023/5_Result/17_scATAC")

####---------Creat ArrowFiles----------
addArchRThreads(threads = 15)
addArchRGenome("hg38")

genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Hsapiens.UCSC.hg38)
genomeAnnotation

inputFiles <- c("/home/liyaru/DATA/MM/MM_ATAC_result/0513/CD19/outs/fragments.tsv.gz")
names(inputFiles) <- c("MM1_BC")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
ArrowFiles

MM1BC <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "MM1_BC",
  copyArrows = TRUE
)

####-------Filter doublets----------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1
)

MM1BC <- filterDoublets(MM1BC)
MM1BC

####------- Filter cells by expression -------------
t <- getMatrixFromProject(
  ArchRProj =MM1BC,
  useMatrix = "GeneScoreMatrix" 
)

tt <- t@assays@data$GeneScoreMatrix

# gene col position
ttt <- rowData(t)
n = which(ttt$name=="CD19")

d <- as.data.frame(t(as.data.frame(tt)))
d[1:5,1:5]

c <- which(d[,n]>0)
CD19_cell <- rownames(d[c,])

MM1BC <- MM1BC[MM1BC$cellNames %in% CD19_cell,]


####------ Run cluster and embedding ------------------
MM1BC <- addIterativeLSI(
  ArchRProj = MM1BC,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force= TRUE
)

MM1BC <- addClusters(
  input = MM1BC,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)

MM1BC <- addUMAP(
  ArchRProj = MM1BC, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

MM1BC <- addImputeWeights(MM1BC)

####-------- Join RNA-------------------
# BC scRNA
BC <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/BC_new.rds")

# scRNA celltype
MM1BC <- addGeneIntegrationMatrix(
  ArchRProj = MM1BC, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  addToArrow = TRUE,
  # use scRNA
  seRNA = BC,
  groupRNA = "celltype",
  # add to scATAC
  nameGroup = "RNA_celltype",
  nameCell = "RNA_celltype_cell",
  nameScore = "RNA_celltype_score",
  force= TRUE
)

# scRNA has_dupli_chr1
MM1BC <- addGeneIntegrationMatrix(
  ArchRProj = MM1BC, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  addToArrow = TRUE,
  # use scRNA
  seRNA = BC,
  groupRNA = "has_dupli_chr1",
  # add to scATAC
  nameGroup = "RNA_dupli_chr1",
  nameCell = "RNA_dupli_chr1_cell",
  nameScore = "RNA_dupli_chr1_score",
  force= TRUE
)


####---- call peak --------------
addArchRThreads(threads = 1)
MM1BC <- addGroupCoverages(ArchRProj = MM1BC, groupBy = "RNA_celltype",force=T)

addArchRThreads(threads = 10)
MM1BC <- addReproduciblePeakSet(
  ArchRProj = MM1BC, 
  groupBy = "RNA_celltype",
  pathToMacs2 = "/home/liyaru/miniconda3/bin/macs2" #PATH
)

MM1BC <- addPeakMatrix(MM1BC)


####---------Differential Peaks ------------------------
markerTest <- getMarkerFeatures(
  ArchRProj = MM1BC, 
  useMatrix = "PeakMatrix",
  groupBy = "RNA_celltype",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Memory-Like 2"
)
# saveRDS(markerTest,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/markerTest.rds")

markerList <- getMarkers(markerTest , cutOff = "FDR < 0.05 & Log2FC > 1")
markerList
diff_peak <- as.data.frame(markerList)
# fwrite(diff_peak,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/diff.peak.csv")

####----Enriched Motif-------------
MM1BC <- addMotifAnnotations(ArchRProj = MM1BC, 
                             motifSet = "cisbp", name = "Motif",force = TRUE)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = MM1BC,
  peakAnnotation = "Motif",
  cutOff = "FDR < 0.05 & Log2FC > 1"
)
motifsUp

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
# fwrite(df,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/diff.peak_enrich.motif.csv")

####--------save and load ArchR proj ----------------------------
pro <- saveArchRProject(ArchRProj = MM1BC)

MM1BC <- loadArchRProject(path = "/media/liyaru/LYR/MM2023/5_Result/17_scATAC/MM1_BC")

####-----------Export peaks-------------------
peak_all <- MM1BC@peakSet
peak_all$RNA_group <- names(peak_all)
names(peak_all) <- paste0(seqnames(peak_all),"_",ranges(peak_all))
peak_all <- as.data.frame(peak_all)
fwrite(peak_all,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/peak_all.csv")

peak_FCRL5 <- peak_all[peak_all$nearestGene =="FCRL5" & !is.na(peak_all$nearestGene),]
fwrite(peak_FCRL5,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/peak_FCRL5.csv")

fwrite(peak_FCRL5[,c("seqnames","start","end")],"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/peak_FCRL5.bed",
       sep = "\t",col.names = F)


####---------open regions-------------------
setwd("/home/liyaru/DATA/MM/MM_ATAC_result/0513/CD19/outs/filtered_peak_bc_matrix")
m <-  Matrix::readMM("matrix.mtx")
b <- fread("barcodes.tsv",header = F)
p <- fread("peaks.bed",header = F)
p$peaks <- paste(p$V1,p$V2,p$V3,sep="_")

mm <- as.matrix(m)
dim(mm)

rownames(mm) <- p$peaks
colnames(mm) <- b$V1
mm[1:5,1:5]

meta <- as.data.frame(MM1BC@cellColData)
meta$Barcode <- rownames(meta)
meta$Barcode2 <- gsub("MM1_BC#","",meta$Barcode)

#meta barcode from ArchR, have a little difference with cellranger peak matrix barcode 
barcode <- b$V1[b$V1 %in% meta$Barcode2]
mmm <- mm[,barcode] %>% as.data.frame()
dim(mmm)

non_zero_counts <- apply(mmm, 2, function(x) sum(x != 0)) %>% as.data.frame()
colnames(non_zero_counts) <- 'Open chromtin regions'
non_zero_counts$cell <- colnames(mmm)

Cell_open_regions <- merge(non_zero_counts,meta,by.x="cell",by.y="Barcode2")

Cell_open_regions 

Cell_open_regions[Cell_open_regions$RNA_celltype == "Memory-Like 2",'Celltype'] <- "Memory-Like 2"
Cell_open_regions[Cell_open_regions$RNA_celltype != "Memory-Like 2",'Celltype'] <- "Other"
fwrite(Cell_open_regions,"/media/liyaru/LYR/MM2023/5_Result/17_scATAC/Open_regions.csv")

p <- ggboxplot(Cell_open_regions,
               x = "Celltype", y = "Open chromtin regions",
               fill = "Celltype")+
  RotatedAxis()+
  scale_fill_manual(values=c("grey","brown3"))+
  NoLegend()
p + stat_compare_means(method = "t.test")

# ggplot(Cell_open_regions, aes(x = "RNA_dupli_chr1", y = "Open chromtin regions")) + 
#   geom_boxplot() +
#   theme_minimal() +
#   labs(title = "Boxplot Example", x = "Groups", y = "Values")
# 
# p <- ggboxplot(Cell_open_regions,
#                x = "RNA_dupli_chr1", y = "Open chromtin regions",
#                fill = "RNA_dupli_chr1")
#   #scale_fill_manual(values=col_CD19)
# p
# 
# Cell_open_regions$RNA_celltype <- factor(Cell_open_regions$RNA_celltype,
#                                          levels = c("Pre B","Pro B","Naive B","Memory-Like 1","Memory-Like 2","Memory-Like 3"))
# p <- ggboxplot(Cell_open_regions,
#                x = "RNA_celltype", y = "Open chromtin regions",
#                fill = "RNA_celltype")+
#   RotatedAxis()
# #scale_fill_manual(values=col_CD19)
# p


###---- Plot QC Figure-----------------
df <- getCellColData(MM1BC, select = c("log10(nFrags)", "TSSEnrichment"))
p1 <- ggPoint(
  x = df[,1],
  y = df[,2],
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

p1

p2 <- plotGroups(
  ArchRProj = MM1BC,
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p2

p3 <- plotGroups(
  ArchRProj = MM1BC,
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)
p3

p4 <- plotFragmentSizes(ArchRProj = MM1BC)
p4

p5 <- plotTSSEnrichment(ArchRProj = MM1BC)
p5


###-------Plot gene score-------------------

p1 <- plotEmbedding(
  ArchRProj = MM1BC, 
  colorBy = "GeneScoreMatrix", 
  name = c("CD19",'CD34','CD38','DNTT','MME'), 
  embedding = "UMAP",
  imputeWeights = NULL,
  plotAs = "points",
  size = 2
)
p1$MME

p2 <- plotEmbedding(
  ArchRProj = MM1BC, 
  colorBy = "GeneScoreMatrix", 
  name = "CD19", 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projHeme3),
  plotAs = "points",
  size = 2
)

p2

###---------Plot embedding----------
p1 <- plotEmbedding(ArchRProj = MM1BC, colorBy = "cellColData", name = "RNA_celltype", embedding = "UMAP",
                    size = 2)
p1

p2 <- plotEmbedding(ArchRProj = MM1BC, colorBy = "cellColData", name = "RNA_dupli_chr1", embedding = "UMAP",
                    size = 2)
p2


p1 + p2

###----------Plot Track-------------
p <- plotBrowserTrack(
  ArchRProj = MM1BC, 
  groupBy = "RNA_Group", 
  geneSymbol = c('CD19','SDC1','CD38','FCRL5'), 
  upstream = 50000,
  downstream = 50000
)
grid::grid.draw(p$CD19)
grid::grid.draw(p$FCRL5)


####------Plot Diff Peak------------
pma <- markerPlot(seMarker = markersPeaks, 
                  name = "Memory-Like 2", 
                  cutOff ="FDR <= 0.05", 
                  plotAs = "MA")
pma

pv <- markerPlot(seMarker = markersPeaks,
                 name = "Memory-Like 2", 
                 cutOff = "FDR <= 0.05", 
                 plotAs = "Volcano")
pv

####--------Plot motif-------------
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp


