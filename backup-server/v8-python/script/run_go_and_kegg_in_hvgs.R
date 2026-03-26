# 加载所需的 R 包
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(patchwork)

# 定义函数：读取交集基因数据并进行 GO/KEGG 富集分析
perform_enrichment_analysis <- function(csv_path, output_dir, is_mouse = TRUE) {
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 读取交集基因的 CSV 文件
  intersection_data <- read.csv(csv_path, stringsAsFactors = FALSE)
  
  # 遍历每个交集列（每列为一个交集基因列表）
  for (column_name in colnames(intersection_data)) {
    # 提取交集基因并去除 NA 和空值
    gene_list <- na.omit(intersection_data[[column_name]])
    gene_list <- gene_list[gene_list != ""]
    
    # 确保基因列表不为空
    if (length(gene_list) == 0) {
      print(paste("Skipping column:", column_name, "- no valid genes found."))
      next
    }
    
    print(paste("Processing column:", column_name, "with", length(gene_list), "genes."))
    
    # GO 富集分析
    go_result <- enrichGO(
      gene          = gene_list,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",       # 选择 GO Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    # 保存 GO 富集结果
    go_output_csv <- file.path(output_dir, paste0(column_name, "_GO_Result.csv"))
    write.csv(go_result@result, go_output_csv, row.names = FALSE)
    print(paste("GO results saved to:", go_output_csv))
    
    # 绘制 GO dotplot
    go_dotplot <- dotplot(go_result, showCategory = 15) +
      ggtitle(paste(column_name, "GO Enrichment"))
    go_dotplot_file <- file.path(output_dir, paste0(column_name, "HVGS_GO_Dotplot.png"))
    go_dotplot_file2 <- file.path(output_dir, paste0(column_name, "HVGS_GO_Dotplot.pDF"))
    ggsave(go_dotplot_file, plot = go_dotplot, height = 8, width = 10)
    ggsave(go_dotplot_file2, plot = go_dotplot, height = 8, width = 10,dpi=500)
    print(paste("GO dotplot saved to:", go_dotplot_file))
    
    
    # KEGG 富集分析
    # 将基因符号转换为 ENTREZ ID
    entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    if (nrow(entrez_ids) == 0) {
      print(paste("Skipping KEGG for column:", column_name, "- no valid ENTREZ IDs found."))
      next
    }
    
    kegg_result <- enrichKEGG(
      gene          = entrez_ids$ENTREZID,
      organism      = 'mmu',     # 小鼠
      keyType       = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05
    )
    
    # 保存 KEGG 富集结果
    kegg_output_csv <- file.path(output_dir, paste0(column_name, "_KEGG_Result.csv"))
    write.csv(kegg_result@result, kegg_output_csv, row.names = FALSE)
    print(paste("KEGG results saved to:", kegg_output_csv))
    
    # 绘制 KEGG dotplot
    kegg_dotplot <- dotplot(kegg_result, showCategory = 15) +
      ggtitle(paste(column_name, "KEGG Enrichment"))
    kegg_dotplot_file <- file.path(output_dir, paste0(column_name, "HVGS_KEGG_Dotplot.png"))
    kegg_dotplot_file2 <- file.path(output_dir, paste0(column_name, "HVGS_KEGG_Dotplot.pdf"))
    ggsave(kegg_dotplot_file, plot = kegg_dotplot, height = 8, width = 10)
    ggsave(kegg_dotplot_file2, plot = kegg_dotplot, height = 8, width = 10,dpi=500)
    print(paste("KEGG dotplot saved to:", kegg_dotplot_file))
  }
}

# 设置输入路径和输出路径
input_csv <- "D:/111/GO_KEGG_Results/HVGS-pairwise_intersections.csv"  # 输入交集基因的 CSV 文件路径
output_directory <- "D:/111/GO_KEGG_Results/intersections"  # 输出文件保存路径

# 运行富集分析
perform_enrichment_analysis(csv_path = input_csv, output_dir = output_directory)
########################################################################################################


library(clusterProfiler)
library(org.Mm.eg.db)  # 如果是人类基因，则改为 org.Hs.eg.db
library(ggplot2)



# 设置输入路径和输出路径
csv_path <- "D:/111/GO_KEGG_Results/unique_hvgs_MG.csv"  # 输入独有基因的 CSV 文件路径
output_dir <- "D:/111/GO_KEGG_Results/unique-MG"  # 输出文件保存路径
# 定义函数：读取独有基因数据并进行 GO/KEGG 富集分析
# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 读取独有基因的 CSV 文件
unique_data <- read.csv(csv_path, stringsAsFactors = FALSE)

# 遍历每个样本
for (sample_name in unique(unique_data$Sample)) {
  # 提取当前样本的基因列表
  gene_list <- unique_data$HVG[unique_data$Sample == sample_name]
  gene_list <- na.omit(gene_list)  # 去除 NA 值
  gene_list <- gene_list[gene_list != ""]  # 去除空值
  
  # 确保基因列表不为空
  if (length(gene_list) == 0) {
    print(paste("Skipping sample:", sample_name, "- no valid genes found."))
    next
  }
  
  print(paste("Processing sample:", sample_name, "with", length(gene_list), "genes."))
  
  # GO 富集分析
  go_result <- enrichGO(
    gene          = gene_list,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",       # 选择 GO Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  # 保存 GO 富集结果
  go_output_csv <- file.path(output_dir, paste0("unique_", sample_name, "_GO_Result.csv"))
  write.csv(go_result@result, go_output_csv, row.names = FALSE)
  print(paste("GO results saved to:", go_output_csv))
  
  # 绘制 GO dotplot
  go_dotplot <- dotplot(go_result, showCategory = 15) +
    ggtitle(paste(sample_name, "GO Enrichment"))
  go_dotplot_file_png <- file.path(output_dir, paste0("unique_", sample_name, "_GO_Dotplot.png"))
  go_dotplot_file_pdf <- file.path(output_dir, paste0("unique_", sample_name, "_GO_Dotplot.pdf"))
  ggsave(go_dotplot_file_png, plot = go_dotplot, height = 8, width = 10)
  ggsave(go_dotplot_file_pdf, plot = go_dotplot, height = 8, width = 10, dpi = 300)
  print(paste("GO dotplot saved to:", go_dotplot_file_png, "and", go_dotplot_file_pdf))
  
  # KEGG 富集分析
  # 将基因符号转换为 ENTREZ ID
  entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  if (nrow(entrez_ids) == 0) {
    print(paste("Skipping KEGG for sample:", sample_name, "- no valid ENTREZ IDs found."))
    next
  }
  
  kegg_result <- enrichKEGG(
    gene          = entrez_ids$ENTREZID,
    organism      = 'mmu',  
    keyType       = "ncbi-geneid",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.3,pvalueCutoff = 0.2
  )
  
  # 保存 KEGG 富集结果
  kegg_output_csv <- file.path(output_dir, paste0("unique_", sample_name, "_KEGG_Result.csv"))
  write.csv(kegg_result@result, kegg_output_csv, row.names = FALSE)
  print(paste("KEGG results saved to:", kegg_output_csv))
  
  # 绘制 KEGG dotplot
  kegg_dotplot <- dotplot(kegg_result, showCategory = 15) +
    ggtitle(paste(sample_name, "KEGG Enrichment"))
  kegg_dotplot_file_png <- file.path(output_dir, paste0("unique_", sample_name, "_KEGG_Dotplot.png"))
  kegg_dotplot_file_pdf <- file.path(output_dir, paste0("unique_", sample_name, "_KEGG_Dotplot.pdf"))
  ggsave(kegg_dotplot_file_png, plot = kegg_dotplot, height = 8, width = 10)
  ggsave(kegg_dotplot_file_pdf, plot = kegg_dotplot, height = 8, width = 10, dpi = 300)
  print(paste("KEGG dotplot saved to:", kegg_dotplot_file_png, "and", kegg_dotplot_file_pdf))
}




