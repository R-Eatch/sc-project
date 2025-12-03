library(Seurat)
library(ggplot2)
library(dplyr)
library(limma)
library(cowplot)
library(wesanderson)
library(scales)

ea <- readRDS("../4.cluster_again/filtered_cca_cleaned.rds")

p1_cca <- DimPlot(ea, reduction = "umap", group.by = "dataset") + ggtitle("CCA - By Dataset")

p2_cca <- DimPlot(ea, reduction = "umap", group.by = "CellType") + ggtitle("CCA - By Cell Type")

p3_cca <- DimPlot(ea, reduction = "umap", group.by = "cluster") + ggtitle("CCA - By Cluster")
ggsave(filename = "dataset.png", plot = p1_cca )
ggsave(filename = "CellType.png", plot = p2_cca,width = 15, height = 10 )
ggsave(filename = "cluster.png", plot = p1_cca )




library(harmony)  
ea <- RunHarmony(ea,group.by.vars="dataset")                              
dataset <- "MG_HARmony" 

ea <- RunUMAP(ea,reduction = "harmony", dims=1:10)

ea <- FindNeighbors(ea, dims = 1:10)
ea <- FindClusters(ea, resolution = 0.1)


p1_cca <- DimPlot(ea, reduction = "umap", group.by = "dataset") + ggtitle("harmony - By Dataset - filtered")

p2_cca <- DimPlot(ea, reduction = "umap", group.by = "CellType") + ggtitle("harmony - By Cell Type - filtered")

ggsave(filename = "p1_harmony.png", plot = p1_cca )
ggsave(filename = "p2_harmony.png", plot = p2_cca,width = 16, height = 10 )

