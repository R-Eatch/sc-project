library(ggplot2)
library(pheatmap)
library(Seurat)
library(reshape2)
library(gridExtra)
library(cowplot)
library(monocle)

# 读取 Seurat 对象
seurat_obj <- readRDS("../1.subset/R-MAGsubset_for_monocle.rds")

# 获取 Seurat 对象中的变量特征基因
variable_genes <- VariableFeatures(seurat_obj)

# 读取 Monocle 对象
cds <- readRDS("../2.monocle/R-MAG_monocle.rds")

# 提取表达矩阵和伪时间信息
expression_matrix <- as.matrix(exprs(cds))
pseudotime <- pData(cds)$Pseudotime
gland <- pData(cds)$gland

# 筛选表达矩阵以仅包含变量特征基因
expression_matrix <- expression_matrix[variable_genes, , drop = FALSE]

# 找到伪时间的最小值和最大值
min_pseudotime <- min(pseudotime)
max_pseudotime <- max(pseudotime)

# 将伪时间分成50个等距段
num_intervals <- 50
pseudotime_breaks <- seq(min_pseudotime, max_pseudotime, length.out = num_intervals + 1)
pseudotime_intervals <- cut(pseudotime, breaks = pseudotime_breaks, include.lowest = TRUE, labels = FALSE)

# 创建一个矩阵来存储每个区间内的相关系数
correlation_matrix <- matrix(NA, nrow = num_intervals, ncol = 1)
rownames(correlation_matrix) <- paste0("Interval_", 1:num_intervals)
colnames(correlation_matrix) <- "Spearman_Correlation"

# 计算每个区间内的 Spearman 相关系数
for (i in 1:num_intervals) {
  # 打印当前进度
  print(paste("Processing interval", i, "of", num_intervals))
  
  # 提取当前区间内的细胞索引
  interval_indices <- which(pseudotime_intervals == i)
  
  # 提取当前区间内的表达矩阵和 gland 信息
  interval_expression_matrix <- expression_matrix[, interval_indices, drop = FALSE]
  interval_gland <- gland[interval_indices]
  
  # 获取当前区间内的两种 gland 的表达矩阵
  gland1_matrix <- interval_expression_matrix[, interval_gland == unique(gland)[1], drop = FALSE]
  gland2_matrix <- interval_expression_matrix[, interval_gland == unique(gland)[2], drop = FALSE]
  
  # 打印调试信息
  print(paste("Interval", i, ": gland1_matrix dimensions:", dim(gland1_matrix)))
  print(paste("Interval", i, ": gland2_matrix dimensions:", dim(gland2_matrix)))
  
  # 检查矩阵是否为空并计算 Spearman 相关系数
  if (!is.null(gland1_matrix) && !is.null(gland2_matrix) && 
      ncol(gland1_matrix) > 1 && ncol(gland2_matrix) > 1) {
    min_cols <- min(ncol(gland1_matrix), ncol(gland2_matrix))
    cor_values <- sapply(1:min_cols, function(k) {
      cor(gland1_matrix[, k], gland2_matrix[, k], method = "spearman", use = "complete.obs")
    })
    correlation_matrix[i, 1] <- mean(cor_values, na.rm = TRUE)
  } else {
    correlation_matrix[i, 1] <- NA
  }
}

# 将相关系数表保存为 CSV 文件
write.csv(correlation_matrix, file = "spearman_correlation_matrix.csv", row.names = TRUE)

correlation_matrix$pseudotime <- 1:50

# 第一张heatmap (基于相关指数绘制)
p1 <- ggplot(data = correlation_matrix, aes(x = pseudotime, y = 1, fill = Spearman_Correlation)) + 
  geom_tile() + 
  scale_fill_gradient(high = "red", low = "white") +
  theme_void() +
  theme(legend.position = "right") +
  ggtitle("Correlation Changes of Rabbit MG and AG Epithelial Cells along Pseudotime")

# 第二张heatmap (基于psetudotime绘制)
p2 <- ggplot(data = correlation_matrix, aes(x = pseudotime, y = 1, fill = pseudotime)) + 
  geom_tile() + 
  scale_fill_gradient(high = "blue", low = "white") +
  theme_void() +
  theme(legend.position = "right") 

# 使用patchwork将两张图组合在一起并添加标题
combined_plot <- p1 / p2 

p3 <- combined_plot

ggsave("spearman_correlation_heatmap_updated.png", plot = p3, width = 15, height = 10)
