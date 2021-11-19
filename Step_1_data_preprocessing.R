###############################################################################################
# Pre-processing for data used to generate Figure 7 and Figure S8 from Kuo, Tasset et al. 2021 

# This script is: STEP 1 of 3

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/Human_PD/data_analysis/DF_2ndRound")

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
T_3922.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_PD/cellranger/T_3922/outs/filtered_feature_bc_matrix")
T_3922 <- CreateSeuratObject(counts = T_3922.counts, project = "Idiopatic_PD_2", min.cells = 3, min.features = 200)
T_3922[["Condition"]] = c('Idiopatic_PD')
T_3922[["Condition_1"]] = c('PD')
T_3922[["Diagnosis"]] = c('Diffuse Lewy Body disease, limbic or transitional type')
T_3922[["Age"]] = c('81')
T_3922[["Sex"]] = c('M')
rm(T_3922.counts)
T_3922[["percent.mt"]] <- PercentageFeatureSet(object = T_3922, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
T_779.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_PD/cellranger/T_779/outs/filtered_feature_bc_matrix")
T_779 <- CreateSeuratObject(counts = T_779.counts, project = "Ctrl_2", min.cells = 3, min.features = 200)
T_779[["Condition"]] = c('Ctrl')
T_779[["Condition_1"]] = c('Ctrl')
T_779[["Diagnosis"]] = c('GBA WT/WT - control')
T_779[["Age"]] = c('76')
T_779[["Sex"]] = c('F')
rm(T_779.counts)
T_779[["percent.mt"]] <- PercentageFeatureSet(object = T_779, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
T_3860.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_PD/cellranger/T_3860/outs/filtered_feature_bc_matrix")
T_3860 <- CreateSeuratObject(counts = T_3860.counts, project = "Familiar_PD_1", min.cells = 3, min.features = 200)
T_3860[["Condition"]] = c('Familiar_PD')
T_3860[["Condition_1"]] = c('PD')
T_3860[["Diagnosis"]] = c('GBA N370S/WT  - heterozygous for GBA')
T_3860[["Age"]] = c('72')
T_3860[["Sex"]] = c('M')
rm(T_3860.counts)
T_3860[["percent.mt"]] <- PercentageFeatureSet(object = T_3860, pattern = "^MT-") #recognize mitochondrial transcripts

#for loading Cell Ranger counts:
T_4432.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/Human_PD/cellranger/T_4432/outs/filtered_feature_bc_matrix")
T_4432 <- CreateSeuratObject(counts = T_4432.counts, project = "Familiar_PD_2", min.cells = 3, min.features = 200)
T_4432[["Condition"]] = c('Familiar_PD')
T_4432[["Condition_1"]] = c('PD')
T_4432[["Diagnosis"]] = c('GBA N370S/WT  - heterozygous for GBA')
T_4432[["Age"]] = c('67')
T_4432[["Sex"]] = c('M')
rm(T_4432.counts)
T_4432[["percent.mt"]] <- PercentageFeatureSet(object = T_4432, pattern = "^MT-") #recognize mitochondrial transcripts

###############################################################################################
all <- T_3922
pdf("T_3922_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("T_3922_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("T_3922_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_3922_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("T_3922_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*13184) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.28, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.28_1107", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("T_3922_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_3922_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.28_1107" #visualizing the singlet vs doublet cells
pdf("T_3922_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"T_3922_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"T_3922_singlets.rds")
singlets<-readRDS("T_3922_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("T_3922_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("T_3922_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"T_3922_singlets_PCA.rds")

###############################################################################################
all <- T_779
pdf("T_779_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("T_779_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("T_779_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_779_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("T_779_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.046*6968) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_321", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("T_779_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_779_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_321" #visualizing the singlet vs doublet cells
pdf("T_779_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"T_779_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"T_779_singlets.rds")
singlets<-readRDS("T_779_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("T_779_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("T_779_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"T_779_singlets_PCA.rds")

###############################################################################################
all <- T_3860
pdf("T_3860_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("T_3860_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("T_3860_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_3860_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("T_3860_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10164) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_772", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("T_3860_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_3860_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.005_772" #visualizing the singlet vs doublet cells
pdf("T_3860_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"T_3860_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"T_3860_singlets.rds")
singlets<-readRDS("T_3860_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("T_3860_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("T_3860_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"T_3860_singlets_PCA.rds")

###############################################################################################
all <- T_4432
pdf("T_4432_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("T_4432_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 & nCount_RNA < 50000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("T_4432_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_4432_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("T_4432_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.061*8872) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.28, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.28, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.28_541", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("T_4432_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("T_4432_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()
Idents(object = all) <- "DF.classifications_0.25_0.28_541" #visualizing the singlet vs doublet cells
pdf("T_4432_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()
saveRDS(all,"T_4432_after_doublet_detection.rds")
#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)
saveRDS(singlets,"T_4432_singlets.rds")
singlets<-readRDS("T_4432_singlets.rds")
Idents(singlets) <- "seurat_clusters"
#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
singlets <- FindVariableFeatures(singlets, selection.method = "vst", nfeatures = 2000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)
#PC capture
pdf("T_4432_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("T_4432_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()
saveRDS(singlets,"T_4432_singlets_PCA.rds")


