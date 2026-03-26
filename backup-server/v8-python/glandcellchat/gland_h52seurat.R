# ==============================================================================
# Step 1: 环境配置
# ==============================================================================
library(reticulate)
library(sceasy)
library(Seurat) # 加载 Seurat 以便进行后续的对象操作

# 【关键】请确保这里的 python 环境里安装了 anndata 和 scanpy
# 如果你知道具体的环境路径，请取消注释并修改下面这一行：
# use_condaenv("/your/path/to/anaconda3/envs/your_env_name", required = TRUE)

# 检查 Python 配置，确保 anndata 可用
py_config()
ad <- import("anndata", convert = FALSE) # 测试是否能加载 anndata

# ==============================================================================
# Step 2: 参数设置 (适配 gland_pp.py)
# ==============================================================================
# 对应 Python 脚本中的 glandlist=['MG','SG','EG']
datasetlist <- c('MG', 'SG', 'EG') 

# 定义输入和输出文件夹
input_dir <- './data'          # gland_pp.py 的输出目录
output_dir <- './data'  # 你的目标输出目录

if(!dir.exists(output_dir)){
  dir.create(output_dir, recursive = TRUE)
}

# ==============================================================================
# Step 3: 转换函数
# ==============================================================================
convert_h5ad_to_seurat <- function(h5ad_path, output_rds_path, project_name) {
  message(paste0("正在转换: ", h5ad_path, " ..."))
  
  tryCatch({
    # 使用 sceasy 将 anndata 转换为 Seurat
    # outFile 参数在某些版本的 sceasy 中可以直接保存，但为了稳妥，我们先转换对象再 saveRDS
    seurat_object <- sceasy::convertFormat(h5ad_path, from="anndata", to="seurat")
     
    # 检查 meta.data 中的列名，确保 Python 中的 celltype/newcelltype 传过来了
    # 如果 Python 里存的是 categorical，这里通常会自动转为 factor
    
    # 保存为 RDS
    saveRDS(seurat_object, file = output_rds_path)
    message(paste0("成功保存: ", output_rds_path))
    
  }, error = function(e) {
    message(paste0("转换失败: ", h5ad_path))
    message(e)
  })
}

# ==============================================================================
# Step 4: 批量执行
# ==============================================================================
for (dataset in datasetlist) {
  # 1. 构建输入路径 (对应 gland_pp.py 的输出: ./data/{dataset}_pp.h5ad)
  h5ad_file <- file.path(input_dir, paste0(dataset, "_pp_subset_for_sceasy.h5ad"))
  
  # 2. 构建输出路径
  rds_file <- file.path(output_dir, paste0(dataset, "_for_cellchat.rds"))
  
  # 3. 执行转换
  if (file.exists(h5ad_file)) {
    convert_h5ad_to_seurat(h5ad_file, rds_file, dataset)
  } else {
    warning(paste0("文件不存在，跳过: ", h5ad_file))
  }
}

message("所有转换任务完成。")
