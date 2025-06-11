library(Seurat)
library(Mfuzz)
library(SingleCellExperiment)
library(dplyr)
cluster <- 6
dataset <- "MRS-MG"      ### M-MG,R-MG,S-MG-,MRS-MG
use_seurat <- TRUE
if (use_seurat){
  
  print("Loading Seurat object...")
  seurat_obj <- readRDS(paste0("D:/111/",dataset,"subset_for_monocle.rds"))
  seurat_obj <- subset(seurat_obj,oldcelltype == 'Luminal cell 1')
  # 从Seurat对象中提取表达数据
  DefaultAssay(seurat_obj) <- 'RNA'
  print("Extracting expression data...")
  expression_data <- GetAssayData(seurat_obj, slot = "data")
  seurat_obj <- FindVariableFeatures(seurat_obj,nfeatures = 4000,selection.method = 'vst')
  expression_data  = expression_data[VariableFeatures(seurat_obj), ]
  # 从Seurat对象中提取时间信息
  print("Extracting time information...")
  time_info <- seurat_obj$time
  
  # 计算每个时间点下每个基因的平均表达水平
  print("Calculating average expression...")
  time_levels <- sort(unique(time_info))
  
  average_expression <- sapply(time_levels, function(t) {
    rowMeans(expression_data[, time_info == t], na.rm = TRUE)
  })
  colnames(average_expression) <- time_levels
  
  
}

# 将平均表达数据转换为MFUZZ格式
print("Converting data to MFUZZ format...")
eset <- ExpressionSet(assayData = as.matrix(average_expression))

# 标准化数据
print("Standardizing data...")
eset <- standardise(eset)

# 处理缺失值
print("Handling missing values...")
exprs(eset)[is.na(exprs(eset))] <- 0

# 选择模糊化参数
print("Estimating fuzzification parameter...")
m <- mestimate(eset)

# 进行模糊C均值聚类
print("Performing fuzzy c-means clustering...")
cl <- mfuzz(eset, c = cluster, m = m)
print("Clustering complete.")

mfuzz.plot2(eset, cl = cl, mfrow = c(2, 5), time.labels = colnames(average_expression))
df <- as.data.frame(cl$cluster)
df1 <- data.frame("gene" = rownames(df),
                  "cluster" = df$`cl$cluster`)
df2 <- as.data.frame(eset@assayData[["exprs"]])
df2$gene <- rownames(df2)
df3 <- merge(df1,df2,by = "gene")
write.csv(df3,paste0(dataset,"Time_gene.csv"),row.names = FALSE)

