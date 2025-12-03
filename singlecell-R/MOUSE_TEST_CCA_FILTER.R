
library(Seurat)
library(ggplot2)
library(dplyr)
library(limma)
library(cowplot)
library(wesanderson)



seurat_combined_cca <- readRDS("MG_CCA_combined.rds")

seurat_filtered <- subset(seurat_combined_cca, subset = CellType != "not_exactly")
p2_cca <- DimPlot(seurat_filtered, reduction = "umap", group.by = "CellType") + ggtitle("CCA - By Cell Type")
ggsave(filename = "1.png", plot = p2_cca, width = 20, height = 10)

dataset <- "filtered_cca"

PC <- 10
res <- 0.2
pK <- 0.05
all_markers <- data.frame()
ea <- seurat_filtered

ea <- NormalizeData(ea)
ea <- FindVariableFeatures(ea, selection.method = "vst", nfeatures = 2000)
ea <- ScaleData(ea, features = rownames(ea))
ea <- RunPCA(ea, features = VariableFeatures(object = ea))
ea <- RunUMAP(ea, dims=1:PC)
ea <- FindNeighbors(ea, dims = 1:PC)
ea <- FindClusters(ea, resolution = res)


p1 <- UMAPPlot(ea, label=T, label.size=8)+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  ggtitle(dataset)
ggsave(plot = p1, filename = "filtered_umap_julei.png" )


# DEG分析
if (length(levels(ea@active.ident)) > 1) {
  markers <- FindAllMarkers(ea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  top5 <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  palette <- wes_palette("Zissou1", 100, type = "continuous")
  p7 <- DoHeatmap(ea, features = as.character(top5$gene), label = TRUE) +
    scale_fill_gradientn(colors = palette)
  
  markers$dataset <- dataset  # 添加样本名作为新列
  
  # 筛选每个cluster的top 10表达基因
  top10_markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  # 将当前样本的top10结果添加到总结果数据框中
  all_markers <- rbind(all_markers, top10_markers)
  
  write.csv(markers, paste0(dataset, "_DEG.csv"))
  
  ggsave(paste0(dataset, "_fig2.pdf"), width = 10, height = 8)
  ggsave(paste0(dataset, "_fig2.tiff"), width = 10, height = 8, dpi = 300)
  
  ea@meta.data[,"Dataset"] <- rep(as.character(dataset), length(row.names(ea@meta.data)))
  saveRDS(ea, paste0(dataset, "_cleaned.rds"))