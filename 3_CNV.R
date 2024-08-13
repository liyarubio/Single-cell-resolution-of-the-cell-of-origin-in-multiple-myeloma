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

setwd("/media/liyaru/LYR/MM2023")
BC <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/BC_new.rds")


#### ---------Run inferCNV (HSC) -----------------
CD34 <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/CD34_new.rds")
Idents(CD34) <- CD34$Celltype
CD34 <- subset(CD34,idents=c("HSC"))
i = CD34
DefaultAssay(i) <- "RNA"

anno <- i@meta.data
anno$cellname <- rownames(anno)
anno <- anno[c('cellname','sample')]

ref_group_names="HD"

counts_matrix <- as.data.frame(GetAssayData(i,assay = "RNA",slot = "count"))
counts_matrix <- counts_matrix[,anno$cellname]

gtf = "/home/liyaru/public_Data/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
gtf = plyranges::read_gff(gtf)

anno_path <- "/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/HSC_celltype_group.xls"
write.table(anno,anno_path, sep = "\t", col.names = F,row.names = F, quote = F)

gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(i)) %>%
  as.data.frame() %>%
  dplyr::select(gene_name, seqnames, start, end) %>%
  dplyr::distinct(gene_name, .keep_all=T) 

gene_path <- paste0("/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/HSC_gene_order_file.xls")
write.table(gene.chr, gene_path, col.names =F, row.names =F, sep = "\t", quote =F )


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix, 
                                    annotations_file=anno_path,
                                    delim="\t",
                                    gene_order_file=gene_path,
                                    ref_group_names=ref_group_names)  

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC",
                              noise_filter = 0.3,
                              num_threads = 18, 
                              cluster_by_groups = T, 
                              cluster_references = F,
                              HMM = T,
                              analysis_mode = "sample",
)


#### ---------Run inferCNV (Plasma) -----------------
MM <- readRDS("/media/liyaru/LYR/MM2023/5_Result/7_Seurat/MM_new.rds")
m = MM@meta.data
m = m[m$cell =='CD138' & m$celltype_new =='Plasmablast',]
plasma <- subset(MM,cells=rownames(m))

i = plasma
DefaultAssay(i) <- "RNA"
table(i$sample)

anno <- i@meta.data
anno$cellname <- rownames(anno)
anno <- anno[c('cellname','sample')]

ref_group_names="HD"

counts_matrix <- as.data.frame(GetAssayData(i,assay = "RNA",slot = "count"))
counts_matrix <- counts_matrix[,anno$cellname]

gtf = "/home/liyaru/public_Data/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
gtf = plyranges::read_gff(gtf)

anno_path <- "/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/P_celltype_group.xls"
write.table(anno,anno_path, sep = "\t", col.names = F,row.names = F, quote = F)

gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(i)) %>%
  as.data.frame() %>%
  dplyr::select(gene_name, seqnames, start, end) %>%
  dplyr::distinct(gene_name, .keep_all=T) 

gene_path <- paste0("/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/P_gene_order_file.xls")
write.table(gene.chr, gene_path, col.names =F, row.names =F, sep = "\t", quote =F )


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix, 
                                    annotations_file=anno_path,
                                    delim="\t",
                                    gene_order_file=gene_path,
                                    ref_group_names=ref_group_names)  

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              out_dir = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/P",
                              noise_filter = 0.3,
                              num_threads = 18, 
                              cluster_by_groups = T, 
                              cluster_references = F,
                              HMM = T,
                              analysis_mode = "sample",
)



#### ---------Run inferCNV (B cells) -----------------
i=BC # B cells seurat object
# same in HSPC and plasma cells

DefaultAssay(i) <- "RNA"
table(i$sample)

anno <- i@meta.data
anno$cellname <- rownames(anno)
anno <- anno[c('cellname','sample')]

ref_group_names="HD"

counts_matrix <- as.data.frame(GetAssayData(i,assay = "RNA",slot = "count"))
counts_matrix <- counts_matrix[,anno$cellname]

gtf = "/home/liyaru/public_Data/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
gtf = plyranges::read_gff(gtf)

anno_path <- "/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/celltype_group.xls"
write.table(anno,anno_path, sep = "\t", col.names = F,row.names = F, quote = F)

gene.chr = gtf %>% plyranges::filter(type == "gene" & gene_name %in% rownames(i)) %>%
  as.data.frame() %>%
  dplyr::select(gene_name, seqnames, start, end) %>%
  dplyr::distinct(gene_name, .keep_all=T) 

gene_path <- paste0("/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/gene_order_file.xls")
write.table(gene.chr, gene_path, col.names =F, row.names =F, sep = "\t", quote =F )


infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts_matrix, 
                                    annotations_file=anno_path,
                                    delim="\t",
                                    gene_order_file=gene_path,
                                    ref_group_names=ref_group_names)  

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff = 0.1,
                              noise_filter = 0.3,
                              out_dir = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC",
                              cluster_by_groups = T, 
                              cluster_references = F,
                              HMM = T,
                              analysis_mode = "sample",
                              num_threads = 18)

#### -------- Add inverCNV result -----------------
BC= infercnv::add_to_seurat(infercnv_output_path="/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC",
                            seurat_obj=BC,top_n=10)


####-----------plot CNV-----------------------
# setwd("/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC")
# setwd("/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC")
setwd("/media/liyaru/LYR/MM2023/5_Result/inferCNV/P")

cnv <- readRDS("run.final.infercnv_obj")

# heatmap
infercnv::plot_cnv(cnv,
                   title = "inferCNV",
                   png_res = 600,
                   output_filename = 'inferCNV',
                   cluster_references = FALSE)

cnv <- readRDS("17_HMM_predHMMi6.hmm_mode-samples.infercnv_obj")

# HMM
infercnv::plot_cnv(cnv,title = "HMM",
                   png_res = 600,
                   output_filename = 'HMM',
                   x.center=3,x.range=c(0,6),
                   cluster_references = FALSE)

#### ------HSC CNV score ---------------
infercnv_obj = readRDS("/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gn <- rownames(expr)

# use when creat inferCNV obj
geneFile <- read.table('/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/HSC_gene_order_file.xls')
# geneFile$w = geneFile$V4 - geneFile$V3    #check no + - chains
sub_geneFile <-  geneFile[geneFile$V1%in%gn,]

df <- as.data.frame(expr) %>% t() %>% as.data.frame()
df$cell <- rownames(df)
df$sample <- substring(rownames(df),1,3)
df[which(df$sample_cnv=="HD_"),'sample_cnv']="HD"

# mean (group X gene)
t=aggregate(df[1:length(df)-1],list(df$sample),mean)

gene_by_sample <- as.data.frame(t(t))
colnames(gene_by_sample) <- gene_by_sample[1,]
gene_by_sample <- gene_by_sample [c(-(1:4)),]

gene_by_sample$gene <- rownames(gene_by_sample)

colnames(sub_geneFile) <- c("gene","chr","start","end")
df<- merge(gene_by_sample,sub_geneFile,by="gene")

df$chr_n <- gsub("chr","",df$chr,)
df$chr_n <- factor(df$chr_n,levels = paste0(1:22))
table(df$chr_n)

df_chr1 <- df[which(df$chr_n=="1"),]
df_chr1$start <- as.numeric(df_chr1$start)
df_chr1 <- df_chr1[order(df_chr1$start),]
fwrite(df_chr1,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/HSC_chr1_CNV_score.csv")

df_chr17 <- df[which(df$chr_n=="17"),]
df_chr17$start <- as.numeric(df_chr17$start)
df_chr17 <- df_chr17[order(df_chr17$start),]
fwrite(df_chr17,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/HSC/HSC_chr17_CNV_score.csv")


#### ------ BCs CNV score---------------
cnv <- BC@meta.data
cnv <- cnv[,c('type',"has_dupli_chr1")]
cnv$cell <- rownames(cnv)

infercnv_obj = readRDS("/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gn <- rownames(expr)

# use when creat inferCNV obj
#geneFile <- read.table('/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/gene_order_file.xls')
sub_geneFile <-  geneFile[geneFile$V1%in%gn,]

df <- as.data.frame(expr) %>% t() %>% as.data.frame()
# rownames(df) <- gsub("HD2","HD",rownames(df))
df$sample <- substring(rownames(df),1,3)
df$cell <- rownames(df)

df <- merge(cnv,df,by="cell")

df[which(df$has_dupli_chr1=="TRUE"),'Dupli']="Chr1_dupli"
df[which(df$has_dupli_chr1=="FALSE"),'Dupli']="Other"
df$sample_cnv <- paste0(df$sample,".",df$Dupli)
table(df$sample_cnv)

df[which(df$sample_cnv=="HD_.Other"),'sample_cnv']="HD"

# check group number
table(df$sample_cnv)

# mean (group X gene)
t=aggregate(df[1:length(df)-1],list(df$sample_cnv),mean)

gene_by_sample <- as.data.frame(t(t))
colnames(gene_by_sample) <- gene_by_sample[1,]
gene_by_sample <- gene_by_sample [c(-(1:4)),]

gene_by_sample$gene <- rownames(gene_by_sample)

colnames(sub_geneFile) <- c("gene","chr","start","end")
df<- merge(gene_by_sample,sub_geneFile,by="gene")

df$chr_n <- gsub("chr","",df$chr,)
df$chr_n <- factor(df$chr_n,levels = paste0(1:22))
table(df$chr_n)

df_chr1 <- df[which(df$chr_n=="1"),]
df_chr1$start <- as.numeric(df_chr1$start)
df_chr1 <- df_chr1[order(df_chr1$start),]
fwrite(df_chr1,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/BC_chr1_CNV_score.csv")


df_chr17 <- df[which(df$chr_n=="17"),]
df_chr17$start <- as.numeric(df_chr17$start)
df_chr17 <- df_chr17[order(df_chr17$start),]
fwrite(df_chr17,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/BC/BC_chr17_CNV_score.csv")

#### ------Plasma CNV score ---------------
infercnv_obj = readRDS("/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/run.final.infercnv_obj")
expr <- infercnv_obj@expr.data
gn <- rownames(expr)

# use when creat inferCNV obj
# geneFile <- read.table('/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/P_gene_order_file.xls')
sub_geneFile <-  geneFile[geneFile$V1%in%gn,]

df <- as.data.frame(expr) %>% t() %>% as.data.frame()
df$cell <- rownames(df)
df$sample <- substring(rownames(df),1,3)
df[which(df$sample_cnv=="HD_"),'sample_cnv']="HD"

# mean (group X gene)
t=aggregate(df[1:length(df)-1],list(df$sample),mean)

gene_by_sample <- as.data.frame(t(t))
colnames(gene_by_sample) <- gene_by_sample[1,]
gene_by_sample <- gene_by_sample [c(-(1:4)),]

gene_by_sample$gene <- rownames(gene_by_sample)

colnames(sub_geneFile) <- c("gene","chr","start","end")
df<- merge(gene_by_sample,sub_geneFile,by="gene")

df$chr_n <- gsub("chr","",df$chr,)
df$chr_n <- factor(df$chr_n,levels = paste0(1:22))
table(df$chr_n)

df_chr1 <- df[which(df$chr_n=="1"),]
df_chr1$start <- as.numeric(df_chr1$start)
df_chr1 <- df_chr1[order(df_chr1$start),]
fwrite(df_chr1,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/P_chr1_CNV_score.csv")

df_chr17 <- df[which(df$chr_n=="17"),]
df_chr17$start <- as.numeric(df_chr17$start)
df_chr17 <- df_chr17[order(df_chr17$start),]
fwrite(df_chr17,file = "/media/liyaru/LYR/MM2023/5_Result/inferCNV/P/P_chr17_CNV_score.csv")



