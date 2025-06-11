# 加载必要的包
library(Seurat)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
ea <- readRDS("../S-MAGCCA_combined.rds")
Idents(ea) <- ea$celltype
#var_features <- VariableFeatures(ea)
var_features <- rownames(ea@assays$RNA)
correlation_results <- list()

# 获取所有时间点和细胞类型
time_points <- unique(ea$time)
cell_types <- unique(ea$celltype)

# 遍历每个time点并计算相关性
for (time_point in time_points) {
  # 提取特定时间点的数据
  time_data <- ea@meta.data %>% filter(time == time_point)
  glands <- unique(time_data$gland)
  
  # 初始化一个临时数据框来存储该时间点的相关性结果
  temp_corr <- data.frame(matrix(NA, nrow = length(cell_types), ncol = 1, dimnames = list(cell_types, "correlation")))
  temp_corr$celltype <- cell_types
  if (length(glands) == 2) {
    # 提取两个腺体的数据
    gland1_cells <- Cells(ea)[which(ea$time == time_point & ea$gland == glands[1])]
    gland2_cells <- Cells(ea)[which(ea$time == time_point & ea$gland == glands[2])]
    
    gland1_obj <- subset(ea, cells = gland1_cells)
    gland2_obj <- subset(ea, cells = gland2_cells)
    
    avg_expr_gland1 <- AverageExpression(gland1_obj, features = var_features, assays = "RNA", slot = "data", return.seurat = FALSE)
    avg_expr_gland2 <- AverageExpression(gland2_obj, features = var_features, assays = "RNA", slot = "data", return.seurat = FALSE)
    
    for (cell_type in cell_types) {
      if (cell_type %in% colnames(avg_expr_gland1[["RNA"]]) && cell_type %in% colnames(avg_expr_gland2[["RNA"]])) {
        gland1_expr <- avg_expr_gland1[["RNA"]][var_features, cell_type, drop = FALSE]
        gland2_expr <- avg_expr_gland2[["RNA"]][var_features, cell_type, drop = FALSE]
        
        if (ncol(gland1_expr) > 0 && ncol(gland2_expr) > 0) {
          corr <- cor(as.numeric(gland1_expr), as.numeric(gland2_expr), method = "spearman")
          temp_corr[cell_type, "correlation"] <- corr
        }
      }
    }
  }
  
  # 将结果添加到correlation_results列表中
  temp_corr$time <- time_point
  correlation_results[[length(correlation_results) + 1]] <- temp_corr
}

# 合并所有相关性结果
correlation_results_df <- do.call(rbind, correlation_results)
correlation_results_df

# 保存为CSV文件
write.csv(correlation_results_df, "gland_correlations.csv", row.names = FALSE)
# 计算细胞类型在不同时间点和腺体中的比例
cell_proportions <- ea@meta.data %>%
  group_by(gland, time, celltype) %>%
  summarise(count = n()) %>%
  group_by(gland, time) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  select(time, gland, celltype, proportion)
# 合并细胞比例和相关性数据
merged_data <- cell_proportions %>%
  left_join(correlation_results_df, by = c("celltype" = "celltype", "time" = "time"))

# 确保时间轴为整数间隔
merged_data$time <- as.integer(merged_data$time)

# 绘制dotplot
p <- ggplot(merged_data, aes(x = time, y = celltype, size = proportion, color = correlation)) +
  geom_point() +
  scale_size_continuous(range = c(1, 10)) +
  scale_color_gradient(low = "green", high = "red") +
  scale_x_continuous(breaks = seq(min(merged_data$time, na.rm = TRUE), max(merged_data$time, na.rm = TRUE), by = 1)) +
  facet_wrap(~gland) +
  theme_minimal() +
  labs(title = "Dotplot of Cell Proportions and Correlations",
       x = "Time",
       y = "Cell Type",
       size = "Cell Proportion",
       color = "Correlation")

# 保存图像
ggsave("dotplot_correlations.png", plot = p, width = 10, height = 8,bg='white')





#############绘制共有基因变化###############
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)

ea <- readRDS("../S-MAGCCA_combined.rds")
# 从 Seurat 对象中提取表达矩阵和元数据
expression_data <- as.data.frame(ea@assays[["RNA"]]@data)
metadata <- ea@meta.data

# 假设表达矩阵的列名与元数据的行名匹配
# 将表达矩阵转置，以便列名变成基因名，行名变成细胞名
expression_data <- t(expression_data)
expression_data <- as.data.frame(expression_data)

# 合并表达矩阵和元数据
combined_data <- cbind(metadata, expression_data)

# 筛选物种为 R 的数据
combined_data <- combined_data %>% filter(species == 'S')

# 去除NA值
combined_data <- na.omit(combined_data)

# 设定表达水平阈值
expression_threshold <- 0.3

# 初始化结果数据框
result <- data.frame(time = integer(), common_genes = integer())

# 按时间点进行分组分析
for (time_point in unique(combined_data$time)) {
  data_at_time <- combined_data %>% filter(time == time_point)
  
  if (nrow(data_at_time) == 0) next
  
  # 分别筛选出两种腺体的数据
  mg_data <- data_at_time %>% filter(gland == 'MG')
  ag_data <- data_at_time %>% filter(gland == 'AG')
  
  if (nrow(mg_data) == 0 | nrow(ag_data) == 0) next
  
  # 找出在两种腺体中都表达的基因
  mg_genes <- colnames(mg_data)[colSums(mg_data[, -c(1:ncol(metadata))] > expression_threshold) > 0]
  ag_genes <- colnames(ag_data)[colSums(ag_data[, -c(1:ncol(metadata))] > expression_threshold) > 0]
  
  common_genes <- intersect(mg_genes, ag_genes)
  
  # 将结果加入到结果数据框中
  result <- rbind(result, data.frame(time = time_point, common_genes = length(common_genes)))
  
  # 打印调试信息
  cat("Time:", time_point, "MG genes:", length(mg_genes), "AG genes:", length(ag_genes), "Common genes:", length(common_genes), "\n")
}

# 绘制并美化折线图
p <- ggplot(result, aes(x = time, y = common_genes)) +
  geom_line(size = 1) +
  labs(title = "The Change of Common Gene Expression Levels Between SugarGlider MG and AG Over Time",
       x = "Time",
       y = "Number of Common Genes") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10)
  )

# 打印图像
print(p)
# 保存图像
ggsave("common_genes_over_time.png", plot = p, width = 10, height = 6)
