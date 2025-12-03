library(Seurat)
library(harmony)
library(cowplot)

datasets <- c(
  "M-MG-E13_5",
  "M-MG-3WK-2",
  "M-MG-LA-1",
  "M-MG-8WK-1",
  "M-MG-LA-2",
  "M-MG-3WK-1",
  "M-MG-GES13_5",
  "M-MG-8WK-2",
  "M-MG-E16-5",
  "M-MG-GES16-5",
  "M-MG-P1"
)


seurat_objects <- list()

for (dataset in datasets) {
  seurat_object_path <- paste0('../2.annotation/', dataset, '_annotation.rds')
  seurat_objects[[dataset]] <- readRDS(seurat_object_path)
}
# 使用CCA进行整合
seurat_objects_cca <- lapply(seurat_objects, function(x) {
  # 确保至少进行了NormalizeData和FindVariableFeatures
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", 
                            nfeatures = 2000, verbose = FALSE)
  return(x)
})

seurat_combined_cca <- merge(seurat_objects_cca[[1]], y = seurat_objects_cca[-1], add.cell.ids = datasets, project = "CCA_combined")

# 执行整合
seurat_combined_cca <- FindIntegrationAnchors(object.list = seurat_objects_cca, dims = 1:20)
seurat_combined_cca <- IntegrateData(anchorset = seurat_combined_cca, dims = 1:20)

# 进行缩放和线性降维
seurat_combined_cca <- ScaleData(seurat_combined_cca)
seurat_combined_cca <- RunPCA(seurat_combined_cca, features = VariableFeatures(object = seurat_combined_cca))

# 为CCA整合运行UMAP
seurat_combined_cca <- RunUMAP(seurat_combined_cca, reduction = "pca", dims = 1:20)

# 绘制UMAP图
p1_cca <- DimPlot(seurat_combined_cca, reduction = "umap", group.by = "dataset") + ggtitle("CCA - By Dataset")
p2_cca <- DimPlot(seurat_combined_cca, reduction = "umap", group.by = "CellType") + ggtitle("CCA - By Cell Type")

combined_cca <- plot_grid(p1_cca, p2_cca, ncol = 2)
ggsave(filename = "CCA_combined_UMAP.png", plot = combined_cca, width = 20, height = 10)

# 保存CCA整合后的Seurat对象
saveRDS(seurat_combined_cca, file = "MG_CCA_combined.rds")

rm(seurat_combined_cca,seurat_objects_cca)


# 使用Harmony进行整合
seurat_objects_harmony <- lapply(seurat_objects, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x)
  x <- RunPCA(x, features = VariableFeatures(object = x))
  return(x)
})

seurat_combined_harmony <- merge(seurat_objects_harmony[[1]], y = seurat_objects_harmony[-1], add.cell.ids = datasets, project = "Harmony_combined")

# 执行Harmony整合
seurat_combined_harmony <- RunHarmony(seurat_combined_harmony, group.by.vars = "dataset", dims = 1:20)



# 为Harmony整合运行UMAP
seurat_combined_harmony <- RunUMAP(seurat_combined_harmony, reduction = "harmony", dims = 1:20)

# 绘制UMAP图
p1_harmony <- DimPlot(seurat_combined_harmony, reduction = "umap", group.by = "dataset") + ggtitle("Harmony - By Dataset")
p2_harmony <- DimPlot(seurat_combined_harmony, reduction = "umap", group.by = "CellType") + ggtitle("Harmony - By Cell Type")

combined_harmony <- plot_grid(p1_harmony, p2_harmony, ncol = 2)
ggsave(filename = "Harmony_combined_UMAP.png", plot = combined_harmony, width = 20, height = 10)

# 保存Harmony整合后的Seurat对象
saveRDS(seurat_combined_harmony, file = "MG_Harmony_combined.rds")

rm(seurat_combined_harmony)

