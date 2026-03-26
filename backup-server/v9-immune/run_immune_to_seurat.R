# 加载 reticulate 包
library(reticulate)

# 指定使用的 Conda 环境的 Python 解释器的绝对路径
#use_condaenv("new_sceasy", required = TRUE)
py_config()
# 加载 sceasy 包
library(sceasy)
path2 <- './forcellchat'
######in python : h5ad_precess.py
datasetlist <- c('M-MG','R-MG','S-MG','R-AG','R-CG','S-AG')
if(!dir.exists(path2)){
  dir.create(path2,recursive = TRUE)
}
# 定义转换函数
convert_h5ad_to_seurat <- function(h5ad_file, output_rds_file) {
  # 使用 sceasy 转换格式
  seurat_object <- sceasy::convertFormat(h5ad_file, from="anndata", to="seurat")

  # 保存 Seurat 对象为 RDS 文件
  saveRDS(seurat_object, file = output_rds_file)
}


#convert_h5ad_to_seurat(paste0('../1.subset/',dataset,'_cleaned.h5ad') , paste0(dataset,"_for_monocle.rds"))
for (dataset in datasetlist){
path1 <- paste0('./',dataset,'/',dataset,'_lum-imm_raw.h5ad')
path3 <- paste0(path2,'/',dataset,'_imm_for_cellchat.rds')
  convert_h5ad_to_seurat(path1 ,path3)
}


