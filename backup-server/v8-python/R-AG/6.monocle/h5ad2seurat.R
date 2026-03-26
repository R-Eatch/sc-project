# 加载 reticulate 包
library(reticulate)

# 指定使用的 Conda 环境的 Python 解释器的绝对路径
#use_condaenv("new_sceasy", required = TRUE)
py_config()
# 加载 sceasy 包
library(sceasy)

######in python : h5ad_precess.py
dataset <- 'R-AG'
# 定义转换函数
convert_h5ad_to_seurat <- function(h5ad_file, output_rds_file) {
  # 使用 sceasy 转换格式
  seurat_object <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

  # 保存 Seurat 对象为 RDS 文件
  saveRDS(seurat_object, file = output_rds_file)
}

# 调用函数进行转换
# 请替换 "path_to_your_h5ad_file.h5ad" 和 "path_to_your_seurat_object.rds" 为你的实际文件路径
#convert_h5ad_to_seurat(paste0('../1.subset/',dataset,'_cleaned.h5ad') , paste0(dataset,"_for_monocle.rds"))
convert_h5ad_to_seurat(paste0(dataset,"_raw.h5ad") , paste0(dataset,"_for_monocle.rds"))
