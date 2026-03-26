#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(reticulate)
  library(sceasy)
})

# -------------------------
# Config
# -------------------------
DATASETS <- c("M-MG","S-MG")

IN_TMPL  <- "./data/%s_sub_cleaned.h5ad"
OUT_DIR  <- "./data/"
USE_COUNTS_AS_X <- FALSE  # 建议 TRUE：确保写入 Seurat 的是 counts

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 如果你需要固定 Python 环境，在这里启用（可选）
# use_condaenv("new_sceasy", required = TRUE)
py_config()

convert_h5ad_to_seurat_rds <- function(h5ad_file, output_rds_file, use_counts_as_X = TRUE) {
  if (!file.exists(h5ad_file)) {
    message("[SKIP] missing: ", h5ad_file)
    return(invisible(NULL))
  }

  tmp_h5ad <- h5ad_file

  if (use_counts_as_X) {
    anndata <- reticulate::import("anndata", convert = FALSE)
    ad <- anndata$read_h5ad(h5ad_file)

    # 要求：X = layers['counts']
    # 若没有 counts layer，会报错；这能帮你尽早发现上游是否正确保存了 counts
    counts <- ad$layers$`__getitem__`("counts")
    ad$X <- counts$copy()

    tmp_h5ad <- tempfile(fileext = ".h5ad")
    ad$write(tmp_h5ad)
  }

  seu <- sceasy::convertFormat(tmp_h5ad, from = "anndata", to = "seurat")
  saveRDS(seu, file = output_rds_file)

  if (use_counts_as_X) {
    unlink(tmp_h5ad)
  }

  message("[OK] ", output_rds_file)
}

# -------------------------
# Run batch
# -------------------------
for (ds in DATASETS) {
  in_h5ad <- sprintf(IN_TMPL, ds, ds)
  out_rds <- file.path(OUT_DIR, sprintf("%s_filtered_reumap.rds", ds))
  convert_h5ad_to_seurat_rds(in_h5ad, out_rds, use_counts_as_X = USE_COUNTS_AS_X)
}

message("Done.")

