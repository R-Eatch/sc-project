library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(enrichplot)
library(ggplot2)

# 加载p8对象
p8 <- readRDS("D:/Rfile/p8.rds")

# 提取p8中的基因和cluster信息
df1 <- p8$annotation_row
df3 <- data.frame("gene" = rownames(df1), "cluster" = df1$Cluster)

# 对每个cluster分别进行GO富集分析
clusters <- unique(df3$cluster)
go_results_list <- list()

for (cluster in clusters) {
  genes <- df3$gene[df3$cluster == cluster]
  ego <- enrichGO(gene         = genes,
                  OrgDb        = org.Mm.eg.db,
                  keyType      = 'SYMBOL',
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  qvalueCutoff  = 0.05)
  ego@result$cluster <- cluster  # 添加cluster信息到结果中
  go_results_list[[as.character(cluster)]] <- ego@result
}

# 合并所有cluster的GO富集结果
all_go_results <- do.call(rbind, go_results_list)

# 提取每个Cluster中p值最小的前5个GO术语
go_top5 <- all_go_results %>%
  group_by(cluster) %>%
  slice_min(order_by = p.adjust, n = 5) %>%
  ungroup() %>%
  arrange(cluster, p.adjust)

# 创建一个包含GO富集结果的data frame，用于绘图
go_df <- data.frame(
  Cluster = factor(go_top5$cluster, levels = unique(go_top5$cluster)),
  GO_term = go_top5$Description,
  p_value = go_top5$p.adjust
)

# 创建GO富集结果的dot plot，不使用pvalue衡量dot的尺寸
dot_plot <- ggplot(go_df, aes(x = factor(1), y = reorder(GO_term, p_value), color = Cluster)) +
  geom_point(size = 5) +
  labs(title = "GO Enrichment Analysis",
       x = "",
       y = "GO Terms") +
  theme_minimal() +
  theme(axis.text.y = element_text(hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# 使用pairwise_termsim计算GO术语之间的相似性
go_results <- enrichGO(gene         = df3$gene,
                       OrgDb        = org.Mm.eg.db,
                       keyType      = 'SYMBOL',
                       ont          = "BP",
                       pAdjustMethod = "BH",
                       qvalueCutoff  = 0.05)

# 计算GO术语之间的相似性
go_results_sim <- pairwise_termsim(go_results)

# 使用enrichplot绘制Enrichment Map
enrichMap <- enrichplot::emapplot(go_results_sim, showCategory = 5) + 
  ggtitle("GO Enrichment Map")

# 保存结果
output_file_dot <- "D:/Rfile/go_enrichment_dotplot.png"
ggsave(output_file_dot, dot_plot, width = 10, height = 8,bg = "white")

output_file_map <- "D:/Rfile/go_enrichment_map.png"
ggsave(output_file_map, enrichMap, width = 10, height = 8)
