library(Seurat)
library(ggplot2)
library(patchwork)

# 读取数据
ea_raw <- readRDS("D:/111/S-MGsubset_for_monocle.rds")
processed_seurat_objects <- list()
PCs <- 50
dataset <- "S-MG"

for (name in unique(ea_raw$rename)){
  print(paste0("Current obj :", name))
  
  # 按照rename字段拆分对象
  ea_processed <- subset(ea_raw, subset = rename == name)
  
  # 标准化数据
  ea_processed <- NormalizeData(ea_processed)
  
  # 寻找变异特征
  ea_processed <- FindVariableFeatures(ea_processed, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  
  # 更新meta数据
  ea_processed@meta.data$orig.ident <- rep(paste0(name, "_", ea_processed@meta.data$orig.ident), nrow(ea_processed@meta.data))
  ea_processed@meta.data$merge_dataset <- rep(name, nrow(ea_processed@meta.data))
  
  # 将处理后的对象加入列表
  processed_seurat_objects <- c(processed_seurat_objects, list(ea_processed))
}

# 选择整合特征
features <- SelectIntegrationFeatures(object.list = processed_seurat_objects)
print("Beginning integration process.")

# 寻找整合锚点
anchors <- FindIntegrationAnchors(object.list = processed_seurat_objects, anchor.features = features, dims = 1:PCs)

# 整合数据
seurat_combined <- IntegrateData(anchorset = anchors, dims = 1:PCs)
print("Data integrated.")

print("Scaling data and running PCA.")
# 进行缩放和线性降维
seurat_combined <- ScaleData(seurat_combined)
seurat_combined <- RunPCA(seurat_combined, features = VariableFeatures(object = seurat_combined))
seurat_combined <- FindNeighbors(seurat_combined, dims = 1:20)
seurat_combined <- FindClusters(seurat_combined, resolution = 0.3)

print("Running UMAP.")
# 进行UMAP降维
ea <- RunUMAP(seurat_combined, reduction = "pca", dims = 1:20)
print("UMAP completed.")
saveRDS(ea,"S-MG-combined.rds")
# 绘制UMAP图
p1 <- DimPlot(ea, reduction = "umap", group.by = "rename") +
  ggtitle(paste0(dataset, " - Subset - By Dataset"))

p2 <- DimPlot(ea, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1) +
  ggtitle(paste0(dataset, " - Subset - By CellType"))

p3 <- DimPlot(ea, reduction = "umap", group.by = "cluster", label = TRUE, pt.size = 1) +
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16, face="bold")) +
  ggtitle(paste0(dataset, " - Subset - By Cluster"))

p4 <- DimPlot(ea, reduction = "umap", group.by = "time") +
  scale_color_viridis_d(option = 'D') +
  ggtitle("UMAP Plot colored by Time")

p5 <- DimPlot(ea, reduction = "umap", group.by = "oldcelltype", label = TRUE, pt.size = 1) +
  ggtitle(paste0(dataset, " - Subset - By OldCellType"))

p6 <- DimPlot(ea, reduction = "umap", group.by = "gland", label = TRUE, pt.size = 1) +
  ggtitle(paste0(dataset, " - Subset - By Gland"))

plotc2 <- (p1 | p2 | p3) / (p4 | p5 | p6)

# 保存图形
ggsave(paste0(dataset, "_CCA.png"), plot = plotc2, width = 30, height = 20)
