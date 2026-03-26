

## -------------------- 0. 依赖 --------------------
suppressPackageStartupMessages({
  library(data.table)   # 读写大型表格
  library(Seurat)       # ≥ v5
  library(dplyr)        # 数据框拼接
})

## -------------------- 1. 文件路径 -----------------
mat.path   <- "D:/Data/GSE109816_normal_heart_umi_matrix.csv.gz"
info1.path <- "D:/Data/GSE109816_normal_heart_cell_info.txt.gz"   # ← 改成你的实际文件
info2.path <- "D:/Data/GSE109816_normal_heart_cell_cluster_info.txt.gz"   # ← 改成你的实际文件

## -------------------- 2. 读表达矩阵 ---------------
ct <- fread(mat.path, data.table = FALSE)
rownames(ct) <- ct[, 1]
ct <- ct[, -1]

## 确保列名（细胞条码）唯一
if (any(duplicated(colnames(ct)))) {
  stop("表达矩阵中检测到重复的细胞 ID，请先去重。")
}

## -------------------- 3. 读两个 cell info ----------
info1 <- fread(info1.path, data.table = FALSE)
info2 <- fread(info2.path, data.table = FALSE)

if (!"ID" %in% names(info1) || !"ID" %in% names(info2)) {
  stop("两个 info 文件都必须包含名为 'ID' 的列。")
}

## -------------------- 4. ID 对齐检查 ---------------
id.mat <- colnames(ct)

missing1 <- setdiff(info1$ID, id.mat)
missing2 <- setdiff(info2$ID, id.mat)

if (length(missing1)) {
  warning(sprintf(
    "info1 中有 %d 个 ID 不在表达矩阵里，将被丢弃：\n%s",
    length(missing1), paste(missing1, collapse = ", ")
  ))
  info1 <- info1[info1$ID %in% id.mat, , drop = FALSE]
}

if (length(missing2)) {
  warning(sprintf(
    "info2 中有 %d 个 ID 不在表达矩阵里，将被丢弃：\n%s",
    length(missing2), paste(missing2, collapse = ", ")
  ))
  info2 <- info2[info2$ID %in% id.mat, , drop = FALSE]
}



## -------------------- 5. 合并元数据 ---------------
meta <- data.frame(ID = id.mat, stringsAsFactors = FALSE) %>%
  left_join(info1, by = "ID", suffix = c("", ".info1")) %>%
  left_join(info2, by = "ID", suffix = c("", ".info2"))

cells.keep <- with(
  meta,
  !is.na(CellType) &              # 真 NA
    !CellType %in% c("NA", "N/A", "Unknown", "")   # 字符串形式
)

meta <- meta[cells.keep, ]
ct   <- ct[, meta$ID]  



rownames(meta) <- meta$ID
meta$ID <- NULL  # Seurat 元数据行名就是细胞条码

## -------------------- 6. 创建 Seurat 对象 ----------
ea <- CreateSeuratObject(
  counts   = ct,
  assay    = "RNA",     # 如有需要可改
  meta.data = meta
)

## -------------------- 7. 基本信息 & 保存 ----------
cat("✓ 对象创建完成：", nrow(ea), "genes ×", ncol(ea), "cells\n")
cat("✓ 元数据列：", paste(colnames(ea@meta.data), collapse = ", "), "\n")

saveRDS(ea, file = "Normal_heart.rds")




ea <- readRDS('D:/Data/Rfile/Normal_heart.rds')
## --------------------------------------------------
## 0. 依赖
## --------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)        # v5.0+
  library(patchwork)     # 拼图用
  library(dplyr)         # 数据操作
})

## 假设 ea 已在工作区；否则 readRDS("ea_seurat.rds")
DefaultAssay(ea) <- "RNA"            # 如有别的 assay，请切换

## --------------------------------------------------
## 1. 归一化 (LogNormalize)
## --------------------------------------------------
ea <- NormalizeData(
  ea,
  normalization.method = "LogNormalize",
  scale.factor         = 1e4,
  verbose = FALSE
)

## --------------------------------------------------
## 2. 高变基因 (HVG) 选取
## --------------------------------------------------
ea <- FindVariableFeatures(
  ea,
  selection.method = "vst",          # 也可 "dispersion"
  nfeatures        = 3000,
  verbose = FALSE
)
hvgs <- VariableFeatures(ea)


## --------------------------------------------------
## 3. 线性缩放 (按细胞均值/方差归一)
##    如果有 UMI 总数 / 线粒体比例等协变量，可放进 vars.to.regress
## --------------------------------------------------
ea <- ScaleData(
  ea,
  features = rownames(ea),           # 全基因缩放
  verbose  = FALSE
)

## --------------------------------------------------
## 4. 主成分分析 (PCA)
## --------------------------------------------------
ea <- RunPCA(
  ea,
  features = VariableFeatures(ea),
  npcs     = 50,                     # 先算 50 PC，后续再截断
  verbose  = FALSE
)

## 快速评估：PCA 可解释方差
ElbowPlot(ea, ndims = 50)            # 选拐点，如 20 PC

## --------------------------------------------------
## 5. 非线性降维 (UMAP @ 20 PCs)
## --------------------------------------------------
npcs.use <- 10
ea <- RunHarmony(ea, "Age")
# 拐点后手动设置
ea <- RunUMAP(
  ea,reduction = "harmony",dims = 1:10
)

## 如果你喜欢 t‑SNE：
# ea <- RunTSNE(ea, dims = 1:npcs.use)

## --------------------------------------------------
## 6. 最近邻图 + Leiden 聚类
## --------------------------------------------------
ea <- FindNeighbors(
  ea,
  dims     = 1:npcs.use,
  k.param  = 20,                     # k ~ 10–30; 依数据规模调节
)

## 多分辨率可尝试 0.2–1.2；这里示例 0.4
ea <- FindClusters(
  ea,
  resolution = 0.2,
)

## --------------------------------------------------
## 7. 可视化
## --------------------------------------------------


## 若之前已存有 celltype，可同时展示
DimPlot(ea, group.by = "CellType")

## --------------------------------------------------
## 8. 统计各 cluster 细胞数
## --------------------------------------------------
cluster_counts <- ea@meta.data %>%
  count(seurat_clusters, name = "n_cells") %>%
  arrange(desc(n_cells))
print(cluster_counts)

## --------------------------------------------------
## 9. 找每个 cluster 的 marker (可选)
## --------------------------------------------------
markers <- FindAllMarkers(
  ea,
  only.pos   = TRUE,
  min.pct    = 0.25,
  logfc.threshold = 0.25
)
head(markers)

## --------------------------------------------------
## 10. 保存
## --------------------------------------------------
saveRDS(ea, file = "D:/Data/NormalHeart_cleaned.rds")


#############################################################################################







## -------------------- 0. 依赖 --------------------
suppressPackageStartupMessages({
  library(data.table)   # 读写大型表格
  library(Seurat)       # ≥ v5
  library(dplyr)        # 数据框拼接
})

## -------------------- 1. 文件路径 -----------------
mat.path   <- "D:/Data/GSE121893_human_heart_sc_umi.csv.gz"
info1.path <- "D:/Data/GSE121893_human_heart_sc_info.txt.gz"   # ← 改成你的实际文件
info2.path <- "D:/Data/GSE121893_all_heart_cell_cluster_info.txt.gz"   # ← 改成你的实际文件

## -------------------- 2. 读表达矩阵 ---------------
ct <- fread(mat.path, data.table = FALSE)
rownames(ct) <- ct[, 1]
ct <- ct[, -1]

## 确保列名（细胞条码）唯一
if (any(duplicated(colnames(ct)))) {
  stop("表达矩阵中检测到重复的细胞 ID，请先去重。")
}

## -------------------- 3. 读两个 cell info ----------
info1 <- fread(info1.path, data.table = FALSE)
info2 <- fread(info2.path, data.table = FALSE)

if (!"ID" %in% names(info1) || !"ID" %in% names(info2)) {
  stop("两个 info 文件都必须包含名为 'ID' 的列。")
}

## -------------------- 4. ID 对齐检查 ---------------
id.mat <- colnames(ct)

missing1 <- setdiff(info1$ID, id.mat)
missing2 <- setdiff(info2$ID, id.mat)

if (length(missing1)) {
  warning(sprintf(
    "info1 中有 %d 个 ID 不在表达矩阵里，将被丢弃：\n%s",
    length(missing1), paste(missing1, collapse = ", ")
  ))
  info1 <- info1[info1$ID %in% id.mat, , drop = FALSE]
}

if (length(missing2)) {
  warning(sprintf(
    "info2 中有 %d 个 ID 不在表达矩阵里，将被丢弃：\n%s",
    length(missing2), paste(missing2, collapse = ", ")
  ))
  info2 <- info2[info2$ID %in% id.mat, , drop = FALSE]
}



## -------------------- 5. 合并元数据 ---------------
meta <- data.frame(ID = id.mat, stringsAsFactors = FALSE) %>%
  left_join(info1, by = "ID", suffix = c("", ".info1")) %>%
  left_join(info2, by = "ID", suffix = c("", ".info2"))

cells.keep <- with(
  meta,
  !is.na(CellType) &              # 真 NA
    !CellType %in% c("NA", "N/A", "Unknown", "")   # 字符串形式
)

meta <- meta[cells.keep, ]
ct   <- ct[, meta$ID]  



rownames(meta) <- meta$ID
meta$ID <- NULL  # Seurat 元数据行名就是细胞条码

## -------------------- 6. 创建 Seurat 对象 ----------
ea <- CreateSeuratObject(
  counts   = ct,
  assay    = "RNA",     # 如有需要可改
  meta.data = meta
)

#hf_types <- c("HF_LA_CM", "HF_LA_NCM", "HF_LV_CM", "HF_LV_NCM")

# 进行子集保留
#ea <- subset(ea, subset = Type %in% hf_types)
## -------------------- 7. 基本信息 & 保存 ----------
cat("✓ 对象创建完成：", nrow(ea), "genes ×", ncol(ea), "cells\n")
cat("✓ 元数据列：", paste(colnames(ea@meta.data), collapse = ", "), "\n")

saveRDS(ea, file = "Disease_heart.rds")




ea <- readRDS('D:/Data/Rfile/Disease_heart.rds')
## --------------------------------------------------
## 0. 依赖
## --------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)        # v5.0+
  library(patchwork)     # 拼图用
  library(dplyr)         # 数据操作
})

## 假设 ea 已在工作区；否则 readRDS("ea_seurat.rds")
DefaultAssay(ea) <- "RNA"            # 如有别的 assay，请切换

## --------------------------------------------------
## 1. 归一化 (LogNormalize)
## --------------------------------------------------
ea <- NormalizeData(
  ea,
  normalization.method = "LogNormalize",
  scale.factor         = 1e4,
  verbose = FALSE
)

## --------------------------------------------------
## 2. 高变基因 (HVG) 选取
## --------------------------------------------------
ea <- FindVariableFeatures(
  ea,
  selection.method = "vst",          # 也可 "dispersion"
  nfeatures        = 3000,
  verbose = FALSE
)
hvgs <- VariableFeatures(ea)


## --------------------------------------------------
## 3. 线性缩放 (按细胞均值/方差归一)
##    如果有 UMI 总数 / 线粒体比例等协变量，可放进 vars.to.regress
## --------------------------------------------------
ea <- ScaleData(
  ea,
  features = rownames(ea),           # 全基因缩放
  verbose  = FALSE
)

## --------------------------------------------------
## 4. 主成分分析 (PCA)
## --------------------------------------------------
ea <- RunPCA(
  ea,
  features = VariableFeatures(ea),
  npcs     = 50,                     # 先算 50 PC，后续再截断
  verbose  = FALSE
)

## 快速评估：PCA 可解释方差
ElbowPlot(ea, ndims = 50)            # 选拐点，如 20 PC

## --------------------------------------------------
## 5. 非线性降维 (UMAP @ 20 PCs)
## --------------------------------------------------
npcs.use <- 10                       # 拐点后手动设置
ea <- RunUMAP(
  ea,
  dims  = 1:npcs.use,
  umap.method = "uwot",
  metric = "cosine",                 # 或 "euclidean"
  verbose = FALSE
)

## 如果你喜欢 t‑SNE：
#ea <- RunTSNE(ea, dims = 1:npcs.use)

## --------------------------------------------------
## 6. 最近邻图 + Leiden 聚类
## --------------------------------------------------
ea <- FindNeighbors(
  ea,
  dims     = 1:npcs.use,
  k.param  = 20,                     # k ~ 10–30; 依数据规模调节
)

## 多分辨率可尝试 0.2–1.2；这里示例 0.4
ea <- FindClusters(
  ea,
  resolution = 0.2,
)

## --------------------------------------------------
## 7. 可视化
## --------------------------------------------------


## 若之前已存有 celltype，可同时展示
DimPlot(ea, group.by = "CellType",reduction = 'umap')
DimPlot(ea, group.by = "group",reduction = 'umap')
## --------------------------------------------------
## 8. 统计各 cluster 细胞数
## --------------------------------------------------
cluster_counts <- ea@meta.data %>%
  count(CellType, name = "n_cells") %>%
  arrange(desc(n_cells))
print(cluster_counts)

## --------------------------------------------------
## 9. 找每个 cluster 的 marker (可选)
## --------------------------------------------------
markers <- FindAllMarkers(
  ea,
  only.pos   = TRUE,
  min.pct    = 0.25,
  logfc.threshold = 0.25
)
head(markers)

## --------------------------------------------------
## 10. 保存
## --------------------------------------------------
#saveRDS(ea, file = "D:/Data/NormalHeart_cleaned.rds")

