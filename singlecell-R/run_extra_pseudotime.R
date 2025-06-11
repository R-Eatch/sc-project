library(Seurat)
library(dplyr)

dataset <- 'R-CG'
use_monocle3 <- TRUE
use_monocle2 <-FALSE
if(use_monocle3){
  library(monocle3)
  cds_subset <- load_monocle_objects(paste0("../6.monocle/",dataset,'_subset_cds.rds'))
  map_df <- data.frame(
    cell_id    = names(pseudotime(cds_subset)),
    pseudotime = pseudotime(cds_subset),
    stringsAsFactors = FALSE
  )
  map_df <- subset(map_df, !is.na(pseudotime))
  
  # 写出 CSV
  out_fn <- "cds_pseudotime_mapping.csv"
  write.csv(map_df, file = out_fn, row.names = FALSE)
  cat(sprintf("monocle3:", out_fn, nrow(map_df)))
}

if(use_monocle2){
  library(monocle)
  cds <- readRDS(paste0("../6.monocle/",dataset,'_monocle.rds'))
  if (!exists("cds"))
    stop("请先加载已有的 CellDataSet 对象到变量 cds")
  
  # 检查 pseudotime 字段名称（新旧版本差异）
  ptime_col <- if ("Pseudotime" %in% colnames(pData(cds))) "Pseudotime" else
    if ("State_Pseudotime" %in% colnames(pData(cds))) "State_Pseudotime" else
      stop("cds 中未找到 Pseudotime 列，请确认已运行 orderCells()")
  map_df <- data.frame(
    cell_id    = rownames(pData(cds)),
    pseudotime = pData(cds)[[ptime_col]],
    stringsAsFactors = FALSE
  )
  # 去掉没有 pseudotime 的细胞（一般是 root 未覆盖）——可选
  map_df <- subset(map_df, !is.na(pseudotime))
  # 写出 CSV
  out_fn <- "cds_pseudotime_mapping.csv"
  write.csv(map_df, file = out_fn, row.names = FALSE)
  cat(sprintf("monocle2:", out_fn, nrow(map_df)))
  
}


