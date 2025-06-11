# 加载必要的R包
library(Seurat)

# 假设你的Seurat对象名为`seurat_object`
seurat_object <- readRDS("3species_v5_noSG.rds")  # 请替换为你的文件路径

# 定义要切割的条件列表
conditions <- list(
  "M-MG" = list(species = "M", gland = "MG"),
  "R-MG" = list(species = "R", gland = "MG"),
  "R-AG" = list(species = "R", gland = "AG"),
  "S-MG" = list(species = "S", gland = "MG"),
  "S-AG" = list(species = "S", gland = "AG"),
  "R-CG" = list(species = "R", gland = "CG"),
)

# 对每个条件进行切割和绘图
for (group_name in names(conditions)) {
  condition <- conditions[[group_name]]
  
  # 筛选Seurat对象
  subset_seurat <- subset(seurat_object, subset = (species == condition$species & gland == condition$gland))
   
  # 保存Seurat对象
  saveRDS(subset_seurat, file = paste0("../newsubset/",group_name"/", group_name, ".rds"))
  
  # 打印保存信息
  print(paste("Saved and Seurat object for", group_name, "in directory", group_name))
}
