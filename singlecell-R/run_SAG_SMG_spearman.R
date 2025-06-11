library(Seurat)
library(pheatmap)
ea1 <- readRDS('D:/111/S-MGsubset_for_monocle_CCA.rds')
ea2 <- readRDS('D:/111/S-AGsubset_for_monocle_CCA.rds')

cluster_mg <- paste0('mg',sort(unique(ea1$seurat_clusters)))
cluster_ag <- paste0('ag',sort(unique(ea2$seurat_clusters)))

corr_df <- data.frame(matrix(
  nrow = length(cluster_mg),
  ncol = length(cluster_ag),
  dimnames = list(cluster_mg,cluster_ag)
  
))
rownames(corr_df) <- cluster_mg
a = 1

for (i in cluster_mg){
  clusteridmg <- gsub('mg','',i)
  ea_subset1 <- subset(ea1,subset = seurat_clusters == clusteridmg )
  df1_aver <- AverageExpression(ea_subset1,assays = 'RNA',layer = 'data',features = VariableFeatures(ea_subset1),verbose = FALSE)
  for (j in cluster_ag){
    clusteridag <- gsub('ag','',j)
    a = a+1
    ea_subset2 <- subset(ea2,subset = seurat_clusters == clusteridag )
    df2_aver <- AverageExpression(ea_subset2,assays = 'RNA',layer = 'data',features = VariableFeatures(ea_subset1),verbose = FALSE)
    corr <- cor(x = as.numeric(df1_aver$RNA),y = as.numeric(df2_aver$RNA),method = 'spearman')
    corr_df[i,j] <- corr
    print(a)
  }
  
}
pheatmap(corr_df,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         display_numbers = TRUE,  # 显示数值
         number_format = "%.2f",  # 保留两位小数
         cluster_rows = TRUE,  # 聚类行
         cluster_cols = TRUE,  # 聚类列
         scale = "none")  # 不进行缩放