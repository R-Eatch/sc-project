library(reticulate)
library(sceasy)
library(Seurat)
library(ggplot2)

# 指定使用的 Conda 环境的 Python 解释器的绝对路径
use_condaenv("new_sceasy", required = TRUE)

# 转换 h5ad 文件为 Seurat 对象
convert_h5ad_to_seurat <- function(h5ad_file, output_rds_file) {
  print(paste("Converting", h5ad_file, "to Seurat object..."))
  seurat_object <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")
  saveRDS(seurat_object, file = output_rds_file)
  print(paste("Saved Seurat object to", output_rds_file))
}



# 主函数
main <- function() {
  input_dir <- "../rawsubset"
  output_dir <- "../newsubset"
  plot_dir <- "../umapplots"
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }

  conditions <- c("M-MG", "R-MG", "R-AG", "S-MG", "S-AG", "M-SG", "S-SG", "R-CG")

  for (condition in conditions) {
    h5ad_file <- file.path(input_dir, paste0(condition, "_raw.h5ad"))
    output_rds_file <- file.path(output_dir,condition,paste0(condition, ".rds"))
    
    convert_h5ad_to_seurat(h5ad_file, output_rds_file)
  }
}

# 运行主函数
main()

