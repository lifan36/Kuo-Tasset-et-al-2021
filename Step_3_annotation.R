###############################################################################################
# Pre-processing for data used to generate Figure 7 and Figure S8 from Kuo, Tasset et al. 2021 

# This script is: STEP 3 of 3

# Adapted from https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html
# by Li Fan 
###############################################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/scratch/lif4001/Human_PD/Inma_4samples")
Human_PD_integrated <- readRDS("Human_PD_integrated_PCA_0.1.rds")
#remove cluster 11, doublets
Human_PD_integrated <- subset(Human_PD_integrated, idents="11", invert=TRUE)

setwd("/athena/ganlab/scratch/lif4001/Human_PD/Inma_4samples/figures")
DefaultAssay(Human_PD_integrated) <- 'RNA'
Human_PD_integrated <- RenameIdents(Human_PD_integrated,
                               `0` = "oligodendrocytes", `1`="astrocytes", `2`="excitatory neurons", `3`="microglia",
                               `4`="inhibitory neurons", `5`="OPCs", `6`="inhibitory neurons", `7`="excitatory neurons",
                               `8`="inhibitory neurons", `9`="inhibitory neurons", `10`="excitatory neurons", `11`="excitatory neurons",
                               `12`="unknown", `13`="oligodendrocytes", `14`="excitatory neurons", `15`="endothelial cells",
                               `16`="excitatory neurons"
)

Human_PD_integrated$celltype <- Idents(Human_PD_integrated)

Idents(Human_PD_integrated) <- "celltype"

pdf("Inma_DotPlot_annotation.pdf", width=9, height=4)
DotPlot(Human_PD_integrated, features = c("EBF1","IGFBP7","VCAN","MYT1","GAD2", "GAD1","CD74","CSF1R","C3","SLC17A7", "CAMK2A", "NRGN","AQP4","GFAP","PLP1", "MBP", "MOBP")) + RotatedAxis()
dev.off()

setwd("/athena/ganlab/scratch/lif4001/Human_PD/Inma_4samples/DEGs")
#Subset celltypes
Cluster_EN <- subset(Human_PD_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(Human_PD_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(Human_PD_integrated, idents = "microglia")
Cluster_AST <- subset(Human_PD_integrated, idents = "astrocytes")
Cluster_OL <- subset(Human_PD_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(Human_PD_integrated, idents = "OPCs")
Cluster_EC <- subset(Human_PD_integrated, idents = "endothelial cells")

Idents(Cluster_EN) <- 'orig.ident'
all_name <- rownames(Cluster_EN)
Cluster_EN_gene_expressions <- AverageExpression(Cluster_EN, features = all_name, assays = 'RNA')
write.csv(Cluster_EN_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_EN_by_samples.csv')

Idents(Cluster_IN) <- 'orig.ident'
all_name <- rownames(Cluster_IN)
Cluster_IN_gene_expressions <- AverageExpression(Cluster_IN, features = all_name, assays = 'RNA')
write.csv(Cluster_IN_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_IN_by_samples.csv')

Idents(Cluster_MG) <- 'orig.ident'
all_name <- rownames(Cluster_MG)
Cluster_MG_gene_expressions <- AverageExpression(Cluster_MG, features = all_name, assays = 'RNA')
write.csv(Cluster_MG_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_MG_by_samples.csv')

Idents(Cluster_AST) <- 'orig.ident'
all_name <- rownames(Cluster_AST)
Cluster_AST_gene_expressions <- AverageExpression(Cluster_AST, features = all_name, assays = 'RNA')
write.csv(Cluster_AST_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_AST_by_samples.csv')

Idents(Cluster_OL) <- 'orig.ident'
all_name <- rownames(Cluster_OL)
Cluster_OL_gene_expressions <- AverageExpression(Cluster_OL, features = all_name, assays = 'RNA')
write.csv(Cluster_OL_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_OL_by_samples.csv')

Idents(Cluster_OPC) <- 'orig.ident'
all_name <- rownames(Cluster_OPC)
Cluster_OPC_gene_expressions <- AverageExpression(Cluster_OPC, features = all_name, assays = 'RNA')
write.csv(Cluster_OPC_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_OPC_by_samples.csv')

Idents(Cluster_EC) <- 'orig.ident'
all_name <- rownames(Cluster_EC)
Cluster_EC_gene_expressions <- AverageExpression(Cluster_EC, features = all_name, assays = 'RNA')
write.csv(Cluster_EC_gene_expressions$RNA, file = 'Inma_average_expression_of_all_genes_in_Cluster_EC_by_samples.csv')







