###############################################################################################
# Pre-processing for data used to generate Figure 7 and Figure S8 from Kuo, Tasset et al. 2021 

# This script is: STEP 2 of 3

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Human_PD/data_analysis/DF_2ndRound")
T_3922 <- readRDS(file = "T_3922_singlets.rds")
T_779 <- readRDS(file = "T_779_singlets.rds")
T_3860 <- readRDS(file = "T_3860_singlets.rds")
T_4432 <- readRDS(file = "T_4432_singlets.rds")

setwd("/athena/ganlab/scratch/lif4001/Human_PD/Inma_4samples")

FPD <- c(T_3860, T_4432)
anchors_FPD <- FindIntegrationAnchors(object.list = FPD, dims = 1:30)
FPD_integrated <- IntegrateData(anchorset = anchors_FPD, dims = 1:30)
rm(T_3860, T_4432, FPD)

Human_PD <- c(T_3922, T_779, FPD_integrated)
anchors_Human_PD <- FindIntegrationAnchors(object.list = Human_PD, dims = 1:30)
Human_PD_integrated <- IntegrateData(anchorset = anchors_Human_PD, dims = 1:30)
rm(IPD_integrated, Ctrl_integrated, FPD_integrated, Human_PD)

pdf("Human_PD_integrated_QC_1.pdf", width=16, height=4)
Idents(Human_PD_integrated) <- "orig.ident"
VlnPlot(object = Human_PD_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()
pdf("Human_PD_integrated_QC_2.pdf", width=16, height=4)
Idents(Human_PD_integrated) <- "Condition"
VlnPlot(object = Human_PD_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(Human_PD_integrated, file = "Human_PD_integrated.rds")

DefaultAssay(Human_PD_integrated) <- 'integrated'

Human_PD_integrated <- ScaleData(Human_PD_integrated, verbose = FALSE)
Human_PD_integrated <- RunPCA(Human_PD_integrated, features = VariableFeatures(object = Human_PD_integrated), verbose = FALSE)

Human_PD_integrated <- FindNeighbors(Human_PD_integrated, dims = 1:20)
Human_PD_integrated <- FindClusters(Human_PD_integrated, resolution = 0.1)
Human_PD_integrated <- RunUMAP(Human_PD_integrated, dims = 1: 20)

str(Human_PD_integrated)

DefaultAssay(Human_PD_integrated) <- 'RNA'
Human_PD_integrated <- NormalizeData(Human_PD_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Human_PD_integrated <- ScaleData(Human_PD_integrated, features = rownames(Human_PD_integrated))

pdf("Human_PD_integrated_umap.pdf", width=6, height=4)
DimPlot(Human_PD_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Human_PD_integrated_umap_split_individual.pdf", width=10, height=6)
DimPlot(Human_PD_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 3)
dev.off()
pdf("Human_PD_integrated_umap_split_Condition.pdf", width=12, height=4)
DimPlot(Human_PD_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

saveRDS(Human_PD_integrated, file = 'Human_PD_integrated_PCA_0.1.rds')

#Human_PD_integrated <- readRDS("Human_PD_integrated_PCA_0.1.rds")

DefaultAssay(Human_PD_integrated) <- 'RNA'
pdf("Human_PD_integrated_umap_test.pdf", width=8, height=6)
DimPlot(Human_PD_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

pdf("Human_PD_integrated_annotation_combine.pdf", width=12, height=6)
sig_all<-c("SYT1","SNAP25","GRIN1","SLC17A7", "CAMK2A", "NRGN","GAD1", "GAD2","PLP1", "MBP", "MOBP","AQP4","GFAP", 
           "CD74","CSF1R","C3","PDGFRA","VCAN","EBF1","IGFBP7","FLT1","CLDN5")
markers.to.plot <- as.matrix(sig_all)
DotPlot(object = Human_PD_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()


