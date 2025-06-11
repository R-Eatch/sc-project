library(clusterProfiler)
library(org.Mm.eg.db)  # 如果是人类基因，则改为 org.Hs.eg.db
library(ggplot2)
library(stringr)

# 设置输入路径和输出路径
input_dir <- "D:/111/go-kegg-monocle2/genes"  # 输入 CSV 文件所在目录
output_dir <- "D:/111/go-kegg-monocle2/"  # 输出文件保存路径

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 定义样本名列表
sample_list <- c(
  'M-MG', 'R-MG', 'S-MG', 'R-AG', 'S-AG',
  'R-CG'
)  # 可根据需要扩展为多个样本

# 获取目录中所有CSV文件
csv_files <- list.files(input_dir, pattern = "*.csv", full.names = TRUE)

# 遍历每个样本
for (sample_name in sample_list) {
  print(paste("Processing sample:", sample_name))
  
  # 筛选与样本名相关的CSV文件
  sample_files <- csv_files[str_detect(csv_files, sample_name)]
  
  for (csv_file in sample_files) {
    # 读取CSV文件
    data <- read.csv(csv_file, stringsAsFactors = FALSE)
    colnames(data)[1:2] <- c('Gene', 'Cluster') 
    # 确保数据包含所需列
    if (!all(c("Gene", "Cluster") %in% colnames(data))) {
      print(paste("Skipping file:", csv_file, "- required columns not found."))
      next
    }
    
    # 遍历每个cluster
    for (cluster_id in unique(data$Cluster)) {
      #if (cluster_id == 1) {
      #  next
      #}
      cluster_data <- data[data$Cluster == cluster_id, "Gene"]
      cluster_data <- na.omit(cluster_data)  # 去除 NA 值
      cluster_data <- cluster_data[cluster_data != ""]  # 去除空值
      
      # 确保基因列表不为空
      if (length(cluster_data) == 0) {
        print(paste("Skipping cluster:", cluster_id, "in file:", csv_file, "- no valid genes found."))
        next
      }
      
      print(paste("Processing cluster:", cluster_id, "in file:", csv_file, "with", length(cluster_data), "genes."))
      
      # GO 富集分析
      go_result <- enrichGO(
        gene          = cluster_data,
        OrgDb         = org.Mm.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",       # 选择 GO Biological Process
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.1,
        qvalueCutoff  = 0.2
      )
      
      # 生成新的文件名
      base_name <- str_sub(basename(csv_file), 5, -5)  # 去掉前3个字符和结尾的.csv
      go_output_csv <- file.path(output_dir, paste0(sample_name, "_", base_name, "_Cluster", cluster_id, "_GO_Result.csv"))
      write.csv(go_result@result, go_output_csv, row.names = FALSE)
      print(paste("GO results saved to:", go_output_csv))
      
      # 绘制并保存 GO dotplot
      plot_title <- paste(sample_name, "Cluster", cluster_id, "GO Enrichment")
      dotplot_file_png <- file.path(output_dir, paste0(sample_name, "_", base_name, "_Cluster", cluster_id, "_GO_Dotplot.png"))
      dotplot_file_pdf <- file.path(output_dir, paste0(sample_name, "_", base_name, "_Cluster", cluster_id, "_GO_Dotplot.pdf"))
      
      # 如果结果为空，生成空白图
      if (nrow(go_result@result) == 0) {
        empty_plot <- ggplot() +
          ggtitle(paste("No GO Results for", sample_name, "Cluster", cluster_id)) +
          theme_void()
        ggsave(dotplot_file_png, plot = empty_plot, height = 8, width = 10)
        ggsave(dotplot_file_pdf, plot = empty_plot, height = 8, width = 10, dpi = 300)
        print(paste("GO results are empty. Saved empty plot to:", dotplot_file_png, "and", dotplot_file_pdf))
      } else {
        # 绘制实际的 dotplot
        dotplot <- dotplot(go_result, showCategory = 15) +
          ggtitle(plot_title)
        ggsave(dotplot_file_png, plot = dotplot, height = 8, width = 10)
        ggsave(dotplot_file_pdf, plot = dotplot, height = 8, width = 10, dpi = 300)
        print(paste("GO dotplot saved to:", dotplot_file_png, "and", dotplot_file_pdf))
      }
    }
  }
}
