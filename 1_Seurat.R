library(Seurat)
library(patchwork)
library(dplyr)
library(tidyr)
library(export)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(ggsci)
library(clustree)
library(monocle)
library(monocle3)
library(data.table)
library(readxl)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(plyr)
library(SCENIC)
library(ggrepel)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(msigdbr)
library(org.Mm.eg.db)
library(infercnv)

####------RDS------------------
MM_CD34 <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/CD34_new.rds")
BC <- readRDS("/media/liyaru/LYR/MM2023/5_Result/1_RDS/6_Bcell_filter.rds")

DimPlot(MM_CD34,group.by = "celltype")
DimPlot(BC1,group.by = "celltype")


####-----Color------------------
colB=c("#90D5E4","#272E6A","#208A42","#F47D2B","#D51F26","#89288F")


####------Read data from cellranger output --------------
setwd("/home/liyaru/DATA/MM")

paths=NULL
mylist=list()

files=c("HD2_CD19","MM1_CD19_s","MM2_CD19",
        "HD2_CD34","MM1_CD34","MM1_CD34_s","MM2_CD34","MM2_CD34_s",
        "HD2_CD138","MM1_CD138_s","MM2_CD138")

sample = c("HD","MM1","MM2",
           "HD","MM1","MM1","MM2","MM2",
           "HD","MM1","MM2")

Cell = c("CD19","CD19","CD19",
         "CD34","CD34","CD34","CD34","CD34",
         "CD138","CD138","CD138")

for (i in c(1:length(files))){
  file=paste0("./result/1_cell_ranger_out/",files[i],"/outs/filtered_feature_bc_matrix")
  paths[i]=file
  
  data<-Read10X(data.dir = file)
  object <- CreateSeuratObject(counts = data,project = files[i],min.cells = 3, min.features = 200)
  
  object[["sample"]] <- sample[i]
  object[["cell"]] <- Cell[i]
  
  mylist<- c(mylist,object)
}
names(mylist)<-files
paths

####--------Merge samples ------------------
merge<- merge(mylist[[1]],
              y=mylist[2:length(files)],
              add.cell.ids=files,
              project="MM")
merge[['percent.mt']] <- PercentageFeatureSet(merge,pattern = "^MT")

####-------- QC (Filter cells) ------------------
# filter cells by basic
merge <- subset(merge,subset = nFeature_RNA>1000 & nFeature_RNA<4000 & nCount_RNA<40000 &  percent.mt<10)

# filter by marker gene expression
t19<-WhichCells(merge, expression = CD19 > 0 & cell=="CD19",slot = "count")
t34<-WhichCells(merge, expression = CD34 > 0 & cell=="CD34" ,slot = "count")
t138<-WhichCells(merge, expression = SDC1 > 0 & cell=="CD138",slot = "count")
cell_exp <-c(t19,t34,t138)

merge_exp <-subset(merge,cells=cell_exp)
table(merge_exp$sample,merge_exp$cell)

####-------- Integrate ( Batch effect )------------------
Idents(merge_exp) <- merge_exp$sample
HD <- subset(merge_exp,ident="HD")
MM1 <- subset(merge_exp,ident="MM1")
MM2 <- subset(merge_exp,ident="MM2")

mylist <- list(HD,MM1,MM2)
mylist2 <- list()

for (i in mylist){
  i[["percent.mt"]] <- PercentageFeatureSet(i,pattern = "^MT")
  i <- subset(i, subset = nFeature_RNA>200 & nFeature_RNA<4000 & percent.mt<10)
  i <- NormalizeData(i, normalization.method = "LogNormalize", scale.factor = 10000)
  i <- FindVariableFeatures(i, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(i)
  i <- ScaleData(i,features = all.genes)
  mylist2 <- c(mylist2,i)
}

anchors <- FindIntegrationAnchors(object.list = mylist2,dims=1:20)
MM <- IntegrateData(anchorset = anchors, dims = 1:20)

####-----Run seurat pipeline------------------
DefaultAssay(MM) <- "integrated"
MM <- MM %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims =1:20) %>%
  RunTSNE(dims = 1:20) %>%
  FindNeighbors() %>%
  FindClusters(res= 0.5) 

DimPlot(MM,group.by = "cell")

MM$cell <- factor(MM$cell,levels = c("CD34","CD19","CD138"))

####-----Annotation with public datasets------------------
# Reference : a dataset of human bone marrow mononuclear (BMNC) cells from eight individual donors, produced by the Human Cell Atlas.
# Stuart*, Butler* et al, Cell 2019
# https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis
# Cell 2019 Comprehensive Integration of Single-Cell Data

# InstallData("bmcite") # load reference data from SeuratData
# bmcite <- LoadData(ds = "bmcite")
load("/home/liyaru/DATA/MM/seurat_data/bmcite.SeuratData_0.3.0/bmcite.SeuratData/data/bmcite.rda")
bmcite

bmcite <- RunUMAP(bmcite, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                  reduction.key = "wnnUMAP_", return.model = TRUE)

p <- DimPlot(bmcite, group.by = "celltype.l2", reduction = "wnn.umap",label=TRUE,repel = TRUE, label.size = 5)+ 
  NoLegend()
p

bmcite<-ScaleData(bmcite,assay = "RNA")
bmcite <- RunSPCA(bmcite, assay = 'RNA', graph = 'wsnn')
bmcite <- FindNeighbors(
  object = bmcite,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE)

DefaultAssay(merge) <- "RNA"
anchors<- FindTransferAnchors(
  reference = bmcite,
  query = MM,
  k.filter = NA,
  reference.reduction = "spca", 
  reference.neighbors = "spca.annoy.neighbors", 
  dims = 1:50
)

MM <- MapQuery(
  anchorset = anchors, 
  query = MM,
  reference = bmcite, 
  refdata = list(
    celltype = "celltype.l2", 
    predicted_ADT = "ADT"),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

DimPlot(MM,group.by = 'predicted.celltype', label.size = 5,label = T)
DimPlot(MM, reduction = 'ref.umap', group.by = 'predicted.celltype', label.size = 5,label = T)

# merge celltype
MM$celltype = MM$predicted.celltype
MM@meta.data[MM@meta.data$predicted.celltype %in% c("CD14 Mono","CD16 Mono","CD4 Memory","CD56 bright NK","cDC2","GMP","Prog_Mk","Prog_RBC"),'celltype'] = "Other"
celltype_level <- c("HSC","LMPP","Prog_B 1","Prog_B 2","Naive B","Memory B","Plasmablast","Other")
MM@meta.data$celltype <- factor(MM@meta.data$celltype,levels = celltype_level)

DimPlot(MM,group.by = "celltype")

####-------HSPC (CD34+) cluster ---------------------------
# HSPCs ref : Dev cell.2022. Temporal molecular program of human hematopoietic stem and progenitor cells after birth
exp = read.csv("/media/liyaru/LYR/MM2023/1_Public_data/6_CD34/hspc_raw_matrix.txt",row.names = 1,sep = "\t")
meta = read.csv("/media/liyaru/LYR/MM2023/1_Public_data/6_CD34/Cell_metadata.txt",row.names = 1,sep = "\t")

# public data process
CD34_ref = CreateSeuratObject(counts = exp)
CD34_ref <- AddMetaData(object = CD34_ref, metadata = meta)
CD34_ref <- CD34_ref %>% 
  NormalizeData()  %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunTSNE(dims = 1:15) %>%
  RunUMAP(dims =1:15)

Idents(MM) <- MM$cell
MM_CD34 <- subset(MM,ident="CD34")
DefaultAssay(MM_CD34) <- "RNA"

anchors <- FindTransferAnchors(reference = CD34_ref, query = MM_CD34,
                               dims = 1:10, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, 
                            refdata = CD34_ref$RNA_clusters,
                            dims = 1:10)

# add reference cell type results to metadata
MM_CD34 <- AddMetaData(MM_CD34, metadata = predictions)

# merge celltypes
table(MM_CD34$predicted.id)
MM_CD34$Celltype = MM_CD34$predicted.id
MM_CD34@meta.data[MM_CD34@meta.data$Celltype %in% c('B_NK1','B_NK2'),'Celltype'] = 'MLP'
MM_CD34@meta.data[MM_CD34@meta.data$Celltype %in% c('Neu1','Neu2','MD','EBM'),'Celltype'] = 'GMP'
MM_CD34@meta.data[MM_CD34@meta.data$Celltype %in% c('Mk','Ery'),'Celltype'] = 'MEP'
MM_CD34@meta.data[MM_CD34@meta.data$Celltype %in% c('MPP1','MPP2'),'Celltype'] = 'MPP'
table(MM_CD34$Celltype)

MM_CD34@meta.data$Celltype = factor(MM_CD34@meta.data$Celltype,levels = c('HSC','MPP','LMPP','MLP','GMP','MEP'))
table(MM_CD34$Celltype)

DefaultAssay(MM_CD34) <- 'integrated'
MM_CD34 <- MM_CD34 %>%
  RunPCA() %>%
  RunTSNE(dims = 1:10) %>%
  RunUMAP(dims =1:10)
DimPlot(MM_CD34,group.by = "Celltype",cols = col_CD34,pt.size = 1)

CD34_celltype  <- as.character(MM_CD34$Celltype)
names(CD34_celltype) <- rownames(MM_CD34@meta.data)

MM$celltype_new <- as.character(MM$celltype)
MM <- AddMetaData(MM,CD34_celltype,col.name = "CD34_celltype")
DimPlot(MM,group.by = "CD34_celltype")

MM@meta.data[!is.na(MM@meta.data$CD34_celltype),"celltype_new"] = MM@meta.data[!is.na(MM@meta.data$CD34_celltype),"CD34_celltype"]


MM$celltype_new <- factor(MM$celltype_new,levels = c('HSC','MPP','LMPP','MLP','GMP','MEP',"Prog_B 1","Prog_B 2","Naive B","Memory B","Plasmablast","Other"))
DimPlot(MM,group.by = "celltype_new")

saveRDS(MM,"/media/liyaru/LYR/MM2023/5_Result/7_Seurat/MM.rds")

####--------B cells (CD19+) cluster-----------------------------
bc = MM@meta.data[MM@meta.data$cell == "CD19" & !MM@meta.data$celltype %in% c("Plasmablast"),] %>% rownames()
BC = subset(MM,cells = bc)

BC <- BC %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims =1:15) %>%
  FindNeighbors() %>%
  FindClusters(res= 0.3)


# m = CD19@reductions$umap@cell.embeddings %>% as.data.frame()
# m = m[m$UMAP_1 <= 8,]   # filter 10 ourlier cells
# BC <- subset(BC,cells=rownames(m))

# MM <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/MM.rds")
# BC <- readRDS("/media/liyaru/LYR/MM2023/5_Result/1_RDS/6_Bcell_filter.rds")

# BMNC label
celltype_BMNC  <- as.character(MM$celltype)
names(celltype_BMNC ) <- rownames(MM@meta.data)
BC <- AddMetaData(BC,celltype_BMNC,col.name = "celltype_BMNC")
BC$celltype_BMNC <- factor(BC$celltype_BMNC,levels=c("Prog_B 1","Prog_B 2","Naive B","Memory B","Other"))
DimPlot(BC,group.by = "celltype_BMNC",pt.size = 0.8)


# Annotation clusters by BMNC celltype and marker genes expression
n=length(unique(BC@meta.data$integrated_snn_res.0.3))
celltype=data.frame(ClusterID=0:(n-1),celltype='Unknown')
celltype[celltype$ClusterID %in% c(0),2]='Naive B'
celltype[celltype$ClusterID %in% c(4),2]='Pro B'
celltype[celltype$ClusterID %in% c(2),2]='Pre B'
celltype[celltype$ClusterID %in% c(1),2]='Memory-Like 1'
celltype[celltype$ClusterID %in% c(3),2]='Memory-Like 2'
celltype[celltype$ClusterID %in% c(5),2]='Memory-Like 3'

for(i in 1:nrow(celltype)){
  BC@meta.data[which(BC@meta.data$integrated_snn_res.0.3 == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

BC$celltype <- factor(BC$celltype,levels = c("Pro B","Pre B","Naive B","Memory-Like 1","Memory-Like 2","Memory-Like 3"))

DimPlot(BC,group.by = "celltype")

Bcelltype  <- as.character(BC$celltype)
names(Bcelltype) <- rownames(BC@meta.data)
MM <- AddMetaData(MM,Bcelltype,col.name = "Bcelltype")
DimPlot(MM,group.by = "Bcelltype")

saveRDS(BC,"/media/liyaru/LYR/MM2023/5_Result/7_Seurat/BC_new.rds")
