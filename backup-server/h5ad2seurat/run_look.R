# 加载Seurat包
library(Seurat)
library(ggplot2)
# 假设你的Seurat对象名为`seurat_object`
seurat_object <- readRDS("3species_v5_noSG.rds")

# 查看元数据
metadata_info <- seurat_object@meta.data
print(head(metadata_info))

# 获取元数据的列名
metadata_columns <- colnames(metadata_info)
print(metadata_columns)

celltype <- table(seurat_object$celltype_new)
gland <- table(seurat_object$gland)
time <- table(seurat_object$time)
species <- table(seurat_object$species)
print(celltype)
print(gland)
print(time)
print(species)
# 绘制DimPlot
p1 <- DimPlot(seurat_object, reduction = "umap",group.by = "celltype_new") + ggtitle("UMAP Plot")

ggsave(p1,filename = "555.png")
