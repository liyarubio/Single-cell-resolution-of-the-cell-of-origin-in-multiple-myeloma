library(clusterProfiler)
library(org.Hs.eg.db)

# remotes::install_github("whtns/seuratTools")
# library("seuratTools")


####-------build seurat-----------------
meta = fread("/media/liyaru/LYR/MM2023/0_DATA/1_Smartseq2/1_Matrix_and_info/information.csv") %>% as.data.frame()
rownames(meta) <- meta$orig.ident

table(meta$replicate,meta$BMPB)

mylist=list()

for (i in c("P1P2","P3P4","P5P6","P7P8")){
  #exp matrix
  exp1 = fread(paste0("/media/liyaru/LYR/MM2023/0_DATA/1_Smartseq2/1_Matrix_and_info/",i,".txt")) %>% as.data.frame()
  
  # transfer ensembl to symbol 50000+features to 30000+features
  gene_id <- bitr(exp1$Geneid, fromType = 'ENSEMBL', toType = c("SYMBOL"), OrgDb = 'org.Hs.eg.db')
  exp1 = merge(gene_id,exp1,by.y="Geneid",by.x="ENSEMBL")
  exp1 <- exp1[!duplicated(exp1$SYMBOL),]
  rownames(exp1) <- exp1$SYMBOL
  exp <- exp1[,8:ncol(exp1)]
  
  # cell names 
  cells = data.frame(colnames(exp))
  colnames(cells) = "full"
  cells = separate(cells,col = "full",into = c("s1","cell","s2"),remove = F,sep="_")
  cells = separate(cells,col = "s1",into = c("s11","s12"),sep="[/]")
  cells$s12 <- gsub("SS","",cells$s12)
  cells$cellname <- paste0("L",cells$s12,".",cells$cell)
  colnames(exp) = cells$cellname
  
  t = intersect(colnames(exp),rownames(meta))
  
  object <- CreateSeuratObject(counts = exp[,t],meta.data=meta)
  
  mylist<- c(mylist,object)
  
}

MMsmart<- merge(mylist[[1]],
              y=mylist[2:4],
              project="MM")

MMsmart@meta.data <- separate(MMsmart@meta.data,col="BMPB",into=c("BMPB","Celltype"),sep = "_")
MMsmart@meta.data$class <- substr(MMsmart@meta.data$replicate, start = 1, stop = 2)
MMsmart$Celltype <- factor(MMsmart$Celltype,levels = c("Pro","Pre","Immature","Breg","NaiveB","Memory B","Plasma"))

MMsmart <- subset(MMsmart, subset = nFeature_RNA>500 & nFeature_RNA<6000 & nCount_RNA < 500000)

VlnPlot(MMsmart,features = c("nCount_RNA","nFeature_RNA"),
        group.by = "replicate",
        pt.size = 0)

VlnPlot(MMsmart,features = c("nCount_RNA","nFeature_RNA"),
        group.by = "Celltype",
        pt.size = 0)

MMsmart <- MMsmart %>%
  NormalizeData() %>%
  FindVariableFeatures %>%
  ScaleData() %>%
  RunPCA() 



####------Run moncole3 HD------------
library(monocle3)
set.seed("618")

Idents(MMsmart) <- MMsmart$class
HD <- subset(MMsmart,ident="HD")

expr <- GetAssayData(HD,assay = "RNA",slot = "count")
# expr <- expr[Bdiff_gene,]

meta = HD@meta.data
gene_anno = data.frame(gene_short_name = rownames(expr))
rownames(gene_anno) = rownames(expr)

cds <- new_cell_data_set(expr,
                         cell_metadata = meta,
                         gene_metadata = gene_anno)

## Normalize and pre-process the data
#cds <- preprocess_cds(cds,num_dim = 82)
cds <- preprocess_cds(cds)

##  Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Cluster the cells
cds <- cluster_cells(cds)

##  Learn a graph
cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by ="Celltype",label_cell_groups = FALSE,
           cell_size = 0.6,group_label_size = 10,label_branch_points = F,
           label_leaves = F,label_roots = F,label_principal_points = F)+
  scale_color_manual(values=col_cds)


# saveRDS(cds,"/media/liyaru/LYR/MM2023/5_Result/1_RDS/monocle_HD.rds")


####------Run moncole3 MM------------
library(monocle3)
set.seed("618")

Idents(MMsmart) <- MMsmart$class
MM <- subset(MMsmart,ident="MM")

expr <- GetAssayData(MM,assay = "RNA",slot = "count")
meta = MM@meta.data
gene_anno = data.frame(gene_short_name = rownames(expr))
rownames(gene_anno) = rownames(expr)

cds <- new_cell_data_set(expr,
                         cell_metadata = meta,
                         gene_metadata = gene_anno)

cds <- preprocess_cds(cds)
#cds <- preprocess_cds(cds,num_dim = 82)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

plot_cells(cds ,color_cells_by ="Celltype",label_cell_groups = FALSE,
           cell_size = 0.6,group_label_size = 10,label_branch_points = F,
           label_leaves = F,label_roots = F,
           label_principal_points = F)+
  scale_color_manual(values=col_cds)

plot_cells(cds ,color_cells_by ="Celltype",label_cell_groups = FALSE,
           cell_size = 0.6,group_label_size = 10,label_branch_points = F,
           label_leaves = F,label_roots = F,
           label_principal_points = F)+
  facet_wrap("~replicate", nrow = 1)+
  scale_color_manual(values=col_cds)

# saveRDS(cds,"/media/liyaru/LYR/MM2023/5_Result/1_RDS/monocle_MM.rds")


####------Run moncole3 all------------
library(monocle3)
set.seed("618")

GOgeneID <- get("GO:0030183", org.Hs.egGO2ALLEGS) 
Bdiff_gene <- mget(GOgeneID, org.Hs.egSYMBOL) %>% unlist() %>% unique()
 
expr <- GetAssayData(MMsmart,assay = "RNA",slot = "count") 
expr <- expr[Bdiff_gene,]

meta = MMsmart@meta.data
gene_anno = data.frame(gene_short_name = rownames(expr))
rownames(gene_anno) = rownames(expr)

cds <- new_cell_data_set(expr,
                         cell_metadata = meta,
                         gene_metadata = gene_anno)

cds <- preprocess_cds(cds,num_dim =8)
#cds <- preprocess_cds(cds,num_dim = 82)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

plot_cells(cds ,color_cells_by ="Celltype",label_cell_groups = FALSE,
           cell_size = 0.6,group_label_size = 10,label_branch_points = F,
           label_leaves = F,label_roots = F,
           label_principal_points = F)+
  scale_color_manual(values=col_cds)

plot_cells(cds ,color_cells_by ="Celltype",label_cell_groups = FALSE,
           cell_size = 0.6,group_label_size = 10,label_branch_points = F,
           label_leaves = F,label_roots = F,
           label_principal_points = F)+
  facet_wrap("~class", nrow = 1)+
  scale_color_manual(values=col_cds)



####------projection ------------
HD <- subset(MMsmart,ident="HD")
MM <- subset(MMsmart,ident="MM")
                       
DefaultAssay(HD) <- "RNA"
DefaultAssay(MM) <- "RNA"

MM <- MM %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()%>%
  RunUMAP(dims = 1:10)

DimPlot(MM,group.by = "Celltype")

HD <- HD %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:10)
DimPlot(HD,group.by = "Celltype")

HD <- RunUMAP(HD,dims = 1:25,return.model=TRUE)

# use monocle UMAP embedding
t = HD@reductions$umap@cell.embeddings
cds.embed <- cds@int_colData$reducedDims$UMAP
t[,1] = cds.embed[rownames(t),1]
t[,2] = cds.embed[rownames(t),2]

HD@reductions$umap@cell.embeddings = t
HD@reductions$umap@misc$model$embedding = t

DimPlot(HD,group.by = "Celltype")

anchors<- FindTransferAnchors(
  reference = HD,
  query = MM,
  k.filter = NA,
  reference.reduction = "pca",
  dims = 1:10
)

MM_project <- MapQuery(
  anchorset = anchors, 
  query = MM,
  reference = HD, 
  refdata = list(celltype = "Celltype"),
  #transferdata.args = list(k.weight = FALSE),
  reference.reduction = "pca", 
  reduction.model = "umap"
)

DimPlot(HD,reduction = "umap",group.by = "Celltype")
DimPlot(MM_project,reduction = "ref.umap",group.by = "Celltype")
#DimPlot(MM_project ,reduction = "ref.umap",group.by = "predicted.celltype")

saveRDS(HD,"/media/liyaru/LYR/MM2023/5_Result/1_RDS/MMsmart_HD.rds")
saveRDS(MM_project,"/media/liyaru/LYR/MM2023/5_Result/1_RDS/MMsmart_MM_project.rds")


# merge HD and MMproject embedding
MMsmart <- RunUMAP(MMsmart,dims = 1:10)
t = MMsmart@reductions$umap@cell.embeddings

he <- HD@reductions$umap@cell.embeddings
me <- MM_project@reductions$ref.umap@cell.embeddings
mm <- rbind(he,me)

t[,1] = mm[rownames(t),1]
t[,2] = mm[rownames(t),2]

MMsmart@reductions$umap@cell.embeddings = t
DimPlot(MMsmart)
saveRDS(MMsmart,"/media/liyaru/LYR/MM2023/5_Result/1_RDS/MMsmart.rds")
