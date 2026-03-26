library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(dplyr)
library(limma)
library(wesanderson)
library(scales)
library(homologene)

########### 定义全局变量 ###############
dataset <- "S-MAG"
datasets <- c("S-MG", "S-AG")  # 下划线替换为连字符
palette <- wes_palette("Darjeeling1", length(datasets), type = "continuous")
use_subset <- FALSE  # 不使用 subset
PCs <- 20

# 读取原始 Seurat 对象
ea1_raw <- readRDS("../S-MG/6.monocle/S-MG_for_monocle.rds")
ea2_raw <- readRDS("../S-AG/6.monocle/S-AG_for_monocle.rds")

# 原始 Seurat 对象列表
ea_raw_list <- list(
  ea1_raw,
  ea2_raw
)

# 初始化一个空列表，用于存储处理后的 Seurat 对象
processed_seurat_objects <- list()

##########################################

# 遍历每个 Seurat 对象，不进行细胞类型筛选
for (i in 1:length(ea_raw_list)) {
  ea_raw <- ea_raw_list[[i]]
  dataset_name <- datasets[i]
  print(paste("Processing dataset:", dataset_name))  # 打印当前正在处理的数据集编号  
  # 直接使用原始对象，不进行 subset
  expression_matrix <- GetAssayData(ea_raw, slot = 'counts')
  metadata <- ea_raw@meta.data
  ea_processed <- CreateSeuratObject(counts = expression_matrix, meta.data = metadata)
  print(paste("Created Seurat object for dataset:", dataset_name))
  processed_seurat_objects[[i]] <- ea_processed
}

for (i in 1:length(processed_seurat_objects)) {
  # 获取当前处理的 Seurat 对象
  current_obj <- processed_seurat_objects[[i]]
  # 不进行时间点筛选，保持注释状态
  # current_obj <- subset(current_obj, timel %in% c(3,4)) 
  # 更新当前对象的所有细胞的 orig.ident 字段为对应的数据集标识符
  current_obj@meta.data$orig.ident <- rep(datasets[i], nrow(current_obj@meta.data))
  current_obj@meta.data$merge_dataset <- rep(datasets[i], nrow(current_obj@meta.data))
  processed_seurat_objects[[i]] <- current_obj
}
print("finish process")

# 合并 Seurat 对象
obj <- Reduce(function(x, y) merge(x, y), processed_seurat_objects)
obj <- JoinLayers(obj)
# 按腺体分割 RNA 数据，替换 $species 为 $gland
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$sample)
print("begin PCA")

# 数据标准化和降维
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = PCs)
print("FINISH PCA")

# 使用 CCA 进行整合
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
print("####################################### FINISH CCA ####################################")

# 聚类和 UMAP 降维
obj <- FindNeighbors(obj, dims = 1:10, reduction = "integrated.cca")
obj <- FindClusters(obj, resolution = 0.2)
print("Running UMAP.")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10)
print("UMAP completed.")

# 保存 Seurat 对象
ea <- JoinLayers(obj)
saveRDS(ea, paste0(dataset, "CCA_combined.rds"))
print("Seurat object saved.")

# 绘制 UMAP 图
p1_cca <- DimPlot(ea, reduction = "umap", group.by = "gland") + 
  scale_colour_manual(values = palette) +
  ggtitle(paste0(dataset, " - CCA - By Dataset"))

p2_cca <- DimPlot(ea, reduction = "umap", group.by = "newcelltype",label = TRUE, label.size = 8) + 
  ggtitle(paste0(dataset, " - CCA - By CellType"))
p4_cca <- DimPlot(ea, reduction = "umap", group.by = "newcelltype",label = TRUE, label.size = 8) +
  ggtitle(paste0(dataset, " - CCA - By CellType"))

p3_cca <- UMAPPlot(ea, label = TRUE, label.size = 8) +
  theme(plot.title = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold")) +
  ggtitle(paste0(dataset, " - CCA - By cluster"))

plotc1 <- p1_cca | p2_cca
plotc2 <- p1_cca | p2_cca | p3_cca
ggsave(paste0(dataset, "_CCA_celltype.png"), plot = plotc1, width = 20, height = 10)
ggsave(paste0(dataset, "_CCA_2.png"), plot = plotc2, width = 30, height = 10)
print("UMAP plots generated and saved.")

# 基于细胞类型寻找差异表达基因
Idents(ea) <- ea$celltype
if (length(levels(ea@active.ident)) > 1) {
  markers <- FindAllMarkers(ea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 筛选每个 cluster 的 top 5 表达基因
  top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  features <- unique(as.character(top5_markers$gene))
  
  # 获取 Zissou1 配色方案
  palette <- wes_palette("Zissou1", 100, type = "continuous")
  
  # DotPlot 使用 Zissou1 配色方案
  p_dot <- DotPlot(ea, features = features) + RotatedAxis() +
    scale_colour_gradientn(colors = palette)

  ggsave(paste0(dataset, "celltype_DotPlot.tiff"), plot = p_dot, width = 15, height = 8, dpi = 300)
  
  # DoHeatmap 使用 Zissou1 配色方案
  p_heatmap <- DoHeatmap(ea, features = features, label = TRUE) +
    scale_fill_gradientn(colors = palette)
  
  ggsave(paste0(dataset, "celltype_Heatmap.tiff"), plot = p_heatmap, width = 10, height = 8, dpi = 100, bg = "white")

  # 保存 DEG 结果
  write.csv(markers, paste0(dataset, "celltype_DEG.csv"))
}

# 基于数据集寻找差异表达基因
Idents(ea) <- ea$merge_dataset
if (length(levels(ea@active.ident)) > 1) {
  markers <- FindAllMarkers(ea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 筛选每个 cluster 的 top 5 表达基因
  top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  features <- unique(as.character(top5_markers$gene))
  
  # 获取 Zissou1 配色方案
  palette <- wes_palette("Zissou1", 100, type = "continuous")
  
  # DotPlot 使用 Zissou1 配色方案
  p_dot <- DotPlot(ea, features = features) + RotatedAxis() +
    scale_colour_gradientn(colors = palette)
  
  ggsave(paste0(dataset, "dataset_DotPlot.tiff"), plot = p_dot, width = 15, height = 8, dpi = 300)
  
  # DoHeatmap 使用 Zissou1 配色方案
  p_heatmap <- DoHeatmap(ea, features = features, label = TRUE) +
    scale_fill_gradientn(colors = palette)
  
  ggsave(paste0(dataset, "dataset_Heatmap.tiff"), plot = p_heatmap, width = 10, height = 8, dpi = 100, bg = "white")

  # 保存 DEG 结果
  write.csv(markers, paste0(dataset, "dataset_DEG.csv"))
}
