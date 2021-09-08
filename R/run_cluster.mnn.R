# run_cluster.mnn.R
# run clusering with MNN-based integration
# Author: Tuo Zhang
# Date: 07/19/2020
# Version: 1.0
# NEW: first version
# 

library(scran)
library(Seurat)
library(dplyr)
library(magrittr)
library(MAST)
library(future)
library(batchelor)
library(scater)
library(ggrepel)
library(ggsci)
library(tibble)

workdir <- "project.folder"
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

refdir <- file.path(sourcedir, "reference")

# project name
project <- "shuibing"

# pattern for defining mitochondrial/ribosomal genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"

# load functions
setwd(workdir)
source("my_functions.R")

# set a random seed
set.seed(98)

# set parallelization for Seurat
plan("multiprocess", workers=6)
# For certain functions, each worker needs access to certain global variables.
# If these are larger than the default limit, you will see errors.
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
dissofile <- file.path(sourcedir, "dissociation", "dissociation_related_genes.human.txt")
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)

# ----------------------------------------------------- filter cells ---------------------------------------------------- #

# sample info
sample.info <- data.frame(Name=c("hPSC-AO"))
rownames(sample.info) <- c("AO")

# load raw UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  sid <- sample.info$Name[k]
  pid <- rownames(sample.info)[k]
  raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid, "filtered_feature_bc_matrix"), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix(raw.counts.list)

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# perform cell filtering
# nGene > 500, nGene <= 6000, nUMI > 1000, nUMI <= 30000, percent.mito < 15%
panc.initial %<>% subset(subset=nFeature_RNA > 500 & nFeature_RNA <= 6000 & nCount_RNA > 1000 & nCount_RNA <= 30000 & percent.mito < 15)

# free memory
rm(raw.counts.list)
rm(raw.counts.all)
rm(tmeta)
rm(tdic)
# ----------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- run MNN-based correction ----------------------------------------------- #
# prepare raw UMI counts table from each donor
selected.donors <- c("AO")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, rownames(subset(panc.initial@meta.data, orig.ident == donor))]
}

# create SingleCellExperiment object
sce.list <- list()
for (donor in names(sample.list)){
  sce.list[[donor]] <- SingleCellExperiment(list(counts=as.matrix(sample.list[[donor]])))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, min.size=200, assay.type="counts", method="hclust", min.mean=0.1))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {normalize(object=x)})

# create a seurat object with raw UMI counts
panc <- CreateSeuratObject(counts=as(do.call(cbind, lapply(sce.list, function(x) counts(x))), "dgCMatrix"), 
                           project=project, assay="RNA", min.cells=0, min.features=0,
                           names.field=1, names.delim="_", meta.data=NULL)

# Calculates the mitochondrial/ribosomal genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)

# Add sample condition
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)

# replace normalized data with the scran normalized data
panc[["RNA"]]@data <- as(do.call(cbind, lapply(sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc), c(disso.genes, 
                                                    grep(mito.pattern, rownames(panc), value=T), 
                                                    grep(ribo.pattern, rownames(panc), value=T)))
VariableFeatures(panc) <- head(variable.genes, 3000)

# scaling the data
panc %<>% ScaleData(features=VariableFeatures(panc))

# Perform linear dimensional reduction
panc %<>% RunPCA(features=VariableFeatures(object=panc))

# free memory
rm(preclust.list)
rm(panc.initial)
rm(sce.list)
rm(tmeta)
rm(tdic)
# ----------------------------------------------------------------------------------------------------------------------- #

# -------------------------------------------- run tSNE/UMAP and clustering --------------------------------------------- #
# select top 21 PCs
pcs <- 21

# Run non-linear dimensional reduction (UMAP/tSNE)
panc %<>% RunUMAP(dims=1:pcs, reduction="pca", n.components=3, seed.use=42, n.neighbors=35, n.epochs=1000)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc %<>% FindNeighbors(reduction="pca", dims=1:pcs)
# FindClusters
panc %<>% FindClusters(resolution=seq(0.05,1,by=0.05), verbose=T)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.5")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# Find differentially expressed features (cluster biomarkers)
for (k in c(0:(length(unique(Idents(panc)))-1))){
  print(paste("C",k,sep=""))
  # wilcox test
  myDETest(panc, k,"wilcox",TRUE,infodir,figdir,tassay="RNA")
}
# ----------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- merge clusters ---------------------------------------------------- #
# save current cluster
panc[["base.clust"]] <- Idents(panc)

# merge clusters after manual reviewing them
# C0 + C1                     ==>  C0 (Ciliated like cells)
# C3 + C4                     ==>  C1 (Goblet like cells)
# C6 + C7 + C8                ==>  C2 (Basal cells)
# C2                          ==>  C3 (Proliferating cells 1)
# C5                          ==>  C4 (Proliferating cells 2)
merge.clust <- c(0,0,3,1,1,4,2,2,2)
names(merge.clust) <- 0:8

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:4)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- final.clust

# Find differentially expressed features (cluster biomarkers)
for (k in c(0:4)){
  print(paste("Cluster",k))
  # wilcox test
  myDETest(panc, k,"wilcox",TRUE,refined.infodir,refined.figdir,tassay="RNA",tsuf="merged_clusters")
}

# save seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))

# free memory
rm(merge.clust)
rm(final.clust)
# ----------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------- generate figures ------------------------------------------------ #
# set x/y-axis boundaries
x.upper <- 6
x.lower <- -16.5
y.upper <- 4
y.lower <- -5.5

# set color for different clusters
my.cluster.color.1 <- c('0'='#fb9a99','1'='#a6cee3','2'='#33a02c','3'='#b2df8a','4'='#1f78b4','5'='#fdbf6f','6'='#ff7f00','7'='#cab2d6','8'='#e31a1c')

# set color for MyFeaturePlot
myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

# Figure 1A
# UMAP plot illustrating five cell populations in the hPSC-AOs
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.cluster.color.1, tlabel=FALSE, tsplit=FALSE, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "UMAPPlot.by_Cluster.png", sep="/"), plot=g, width=10, height=8, dpi=300)

# Figure 1B
# UMAP and violin plots showing the expression of selected genes
for (gene in c('FOXJ1','MUC5B','KRT5','MKI67')){
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tncol=1, tlegend=NULL, tptsize=1, talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:4, tcolor_by='ident', tcolor_order=0:4, tcolor=my.cluster.color.1, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"png",sep="."), sep="/"), plot=plot, height=4, width=6, dpi=300)
}

# Figure 1F
# UMAP and violin plots showing the expression of SARS-CoV-2 entry factors
for (gene in c('ACE2','CTSL','FURIN','TMPRSS2','NRP1')){
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tncol=1, tlegend=NULL, tptsize=1, talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:4, tcolor_by='ident', tcolor_order=0:4, tcolor=my.cluster.color.1, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"png",sep="."), sep="/"), plot=plot, height=4, width=6, dpi=300)
}

# Figure S2A & B
# UMAP showing the expression levels of (A) hPSC-AOs cell markers and (B) alveolar epithelial type 2 cell markers
for (gene in c('CDK1','TOP2A','TP63','KRT17','SFTPD','SFTPC','NAPSA')){
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, tncol=1, tlegend=NULL, tptsize=1, talpha=0.7, twidth=15, theight=12.5, tunits="in", tres=300)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"png",sep='.')), plot=plot, width=10, height=8, dpi=300)
}

# load the reference marker genes
gene.list <- read.table(file.path(sourcedir, 'reference', 'ref.genes.txt.gz'), header=T, check.names=F, stringsAsFactors=F, sep='\t')
ref.markers <- list()
for (c in unique(gene.list$Cluster)){
  ref.markers[[paste0('Ref_',c)]] <- as.vector(subset(gene.list, Cluster == c)[,'Gene'])
}

# load the marker genes for our clusters
our.markers <- list()
for (k in 0:8){
  print(paste("load","our","cluster",k))
  # load marker
  tmarker <- read.table(file.path(infodir, paste0("C",k,".wilcox.origin.markers.txt")), header=T, check.names=F, stringsAsFactors=F, row.names=1, sep='\t')
  # filter based on p_val_adj and avg_logFC
  tmarker <- subset(tmarker, p_val_adj < padj.cut & avg_logFC > 0)
  our.markers[[paste0('Our_',k)]] <- rownames(tmarker)
}

# Figure 1D
# Correlation analysis of genes with cell fates in hPSC-AOs and adult human lung cells
plotFracOverlap(ref.markers, our.markers, infodir, figdir, "correlation.ref_cluster", tsvg=F)

# Figure 1C
# Enrichement analysis of hPSC-AOs using genes highly expressed in adult human ciliates or proximal ciliated cells
scoreByRefMarkers(panc, ref.markers, trefPattern='^Ref_', infodir, figdir, 'enrichment.ref_cluster', tsvg=F)
# ----------------------------------------------------------------------------------------------------------------------- #
