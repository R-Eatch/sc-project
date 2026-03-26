#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

# ========= 可配置参数 =========
INPUT_RDS       <- "../data/splits_seurat/combined_seurat_full.rds"
OUTDIR          <- "./result"                  # 输出根目录：result/{group-name}/...
AGE_COL         <- "Age_group"
TYPE_COL        <- "Main_cell_type"
COMBO_COL_NAME  <- "age_celltype"

SOURCE_ASSAY    <- "RNA"       # 从哪个 assay 读取 counts/data 来构建 RNA_sym
ASSAY_USE       <- "RNA_sym"   # 后续分析与绘图使用的 assay

TEST_USE        <- "wilcox"
LOGFC_THRESH    <- 0.25
MIN_PCT         <- 0.10
PADJ_THRESH     <- 0.05
TOPN_PER_GROUP  <- 10

UMAP_WIDTH      <- 9
UMAP_HEIGHT     <- 7
DOT_WIDTH       <- 10
DOT_HEIGHT      <- 8

# ========= 小工具 =========
ensure_dir <- function(d) { if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE) }

save_plot <- function(p, file, w, h) {
  ggsave(file, p, width = w, height = h, dpi = 300, limitsize = FALSE)
  message("Saved: ", file)
}

safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[/\\\\:?*\"<>|]", "_", x, perl = TRUE)  # 非法文件名字符替换
  x <- gsub("\\s+", "_", x, perl = TRUE)             # 连续空白替换为下划线
  ifelse(nchar(x) == 0, "NA", x)
}

# ========= 构建 RNA_sym：用 gene_name 重命名行，并复用已标准化 data 槽 =========
# - 从 SOURCE_ASSAY 的 counts 与 data 读取
# - 用 meta.features$gene_name（缺失则回退行名）生成符号名
# - 去空白/NA，唯一化；再按你的习惯把 '_' 替换为 '-'
# - 新建 RNA_sym，写入 counts 与 data；meta.features 保存 gene_id 与 gene_name
build_RNA_sym <- function(seurat_obj, source_assay = SOURCE_ASSAY, target_assay = "RNA_sym",
                          gene_name_col = "gene_name", replace_underscore_with_dash = TRUE) {
  if (!source_assay %in% Assays(seurat_obj)) {
    stop("Source assay '", source_assay, "' not found. Available: ", paste(Assays(seurat_obj), collapse = ", "))
  }
  ass <- GetAssay(seurat_obj, assay = source_assay)

  # 原始 IDs 与 counts / data
  cts <- GetAssayData(seurat_obj, assay = source_assay, slot = "counts")
  if (is.null(cts) || nrow(cts) == 0) stop("No counts in assay ", source_assay)
  gene_id <- rownames(cts)

  rna_data <- tryCatch(GetAssayData(seurat_obj, assay = source_assay, slot = "data"),
                       error = function(e) NULL)
  if (is.null(rna_data) || nrow(rna_data) != nrow(cts) || ncol(rna_data) != ncol(cts)) {
    stop("No usable normalized 'data' slot to reuse in assay ", source_assay,
         ". Please normalize upstream (e.g., NormalizeData) before building RNA_sym.")
  }

  # 取 gene_name 映射（优先 meta.features$gene_name；否则回退行名）
  mf <- ass@meta.features
  if (!is.null(mf) && nrow(mf) > 0 && gene_name_col %in% colnames(mf)) {
    sym <- mf[gene_id, gene_name_col, drop = TRUE]
  } else {
    sym <- gene_id
  }
  sym <- as.character(sym)
  sym <- trimws(sym)
  fallback_ids <- gene_id
  mask <- is.na(sym) | sym == "" | sym == "NA"
  sym[mask] <- fallback_ids[mask]
  sym <- trimws(sym)
  mask2 <- is.na(sym) | sym == ""
  sym[mask2] <- fallback_ids[mask2]
  stopifnot(length(sym) == length(gene_id))

  # 唯一化；按需把下划线替换为短横（与你给的片段一致）
  sym_unique <- make.unique(sym)
  if (replace_underscore_with_dash) {
    sym_clean <- gsub("_", "-", sym_unique, fixed = TRUE)
  } else {
    sym_clean <- sym_unique
  }

  # 应用到 counts/data 行名
  rownames(cts)      <- sym_clean
  rownames(rna_data) <- sym_clean

  # 构建 target assay，并写入 meta.features（保存原 gene_id 与 gene_name）
  options(Seurat.object.assay.version = "v5")
  new_assay <- CreateAssayObject(counts = cts)
  new_assay@data <- rna_data[rownames(new_assay@counts), colnames(seurat_obj), drop = FALSE]
  mf_new <- data.frame(
    gene_id   = gene_id,
    gene_name = sym_clean,
    row.names = sym_clean,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  new_assay@meta.features <- mf_new

  seurat_obj[[target_assay]] <- new_assay
  DefaultAssay(seurat_obj) <- target_assay
  seurat_obj
}

# ========= feature 映射与作图/导出 =========
# 从指定 assay 抽取 gene_id / gene_name 映射（若缺失则用行名）
get_feature_map <- function(seu, assay = DefaultAssay(seu)) {
  ass <- GetAssay(seu, assay = assay)
  fid <- rownames(ass)
  fm  <- ass@meta.features
  # 若 meta.features 带 gene_id/gene_name（RNA_sym 会有），直接使用
  if (!is.null(fm) && nrow(fm) > 0 &&
      all(c("gene_id", "gene_name") %in% colnames(fm))) {
    tibble(gene_id = fm$gene_id, gene_name = fm$gene_name)
  } else {
    # 兜底：gene_id = 行名；gene_name 也用行名
    tibble(gene_id = fid, gene_name = fid)
  }
}

select_top_markers <- function(markers_df, topn = TOPN_PER_GROUP, padj = PADJ_THRESH) {
  markers_df %>%
    filter(p_val_adj < padj) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE) %>%
    slice_head(n = topn) %>%
    ungroup()
}

# DotPlot：内部以行名（RNA_sym 的行名即 gene_name）取数；x 轴显示 gene_name
# 若你后续把行名设置为别的 ID，这里也可用 feature_map 重映射
dotplot_from_markers <- function(seu, features_df, group_col, assay = ASSAY_USE) {
  # features_df 需含 gene_name（此时 RNA_sym 的行名就是 gene_name）
  ord_names <- unique(features_df$gene_name)
  p <- DotPlot(seu, features = ord_names, group.by = group_col, assay = assay) +
       RotatedAxis()
  p
}

umap_plot_one_group <- function(seu, group_col, out_dir) {
  nlev <- nlevels(factor(seu[[group_col]][,1]))
  add_label <- nlev <= 50
  p <- DimPlot(
    seu, reduction = "umap",
    group.by = group_col,
    label = add_label, repel = add_label, label.size = 3
  )
  save_plot(p, file.path(out_dir, paste0("UMAP_", group_col, ".png")), UMAP_WIDTH, UMAP_HEIGHT)
  save_plot(p, file.path(out_dir, paste0("UMAP_", group_col, ".pdf")), UMAP_WIDTH, UMAP_HEIGHT)
}

# FeaturePlot：RNA_sym 的行名已经是 gene_name，直接传 gene_name 作图；每小面板标题也用 gene_name
featureplot_top_per_identity <- function(seu, top_markers, group_col, out_dir_feature, assay = ASSAY_USE) {
  identities <- unique(top_markers$cluster)
  ensure_dir(out_dir_feature)
  for (idv in identities) {
    df <- top_markers %>% filter(cluster == idv)
    if (nrow(df) == 0) next
    genes_name <- unique(df$gene_name)

    plist <- FeaturePlot(
      seu, features = genes_name, reduction = "umap",
      combine = FALSE, order = TRUE, cols = c("lightgrey", "red")
    )
    for (i in seq_along(plist)) {
      plist[[i]] <- plist[[i]] + ggtitle(genes_name[i])
    }
    w <- min(12, max(8, ceiling(length(genes_name) / 2) * 4))
    h <- min(10, max(6, ceiling(length(genes_name) / 3) * 3))
    p_comb <- wrap_plots(plist, guides = "collect") + plot_annotation(
      title = paste0("FeaturePlot top", length(genes_name), " genes for ", group_col, " = ", idv)
    )

    fn_base <- paste0("FeaturePlot_", group_col, "_", safe_name(idv), "_top", length(genes_name))
    save_plot(p_comb, file.path(out_dir_feature, paste0(fn_base, ".png")), w, h)
    save_plot(p_comb, file.path(out_dir_feature, paste0(fn_base, ".pdf")), w, h)
  }
}

run_deg_one_group <- function(seu, group_col, out_dir_group, feature_map, assay = ASSAY_USE) {
  Idents(seu) <- factor(seu[[group_col]][,1])  # 以该列为身份

  # FindAllMarkers：所有组 vs 其它；这里会在 ASSAY_USE（RNA_sym）上运行
  markers <- FindAllMarkers(
    object = seu, only.pos = TRUE, test.use = TEST_USE,
    logfc.threshold = LOGFC_THRESH, min.pct = MIN_PCT, assay = assay
  )

  # 将 Seurat 的 gene 列（此时即 RNA_sym 的行名 = gene_name）与 feature_map 合并
  markers <- markers %>%
    rename(gene_name = gene) %>%
    left_join(feature_map, by = "gene_name") %>%
    relocate(gene_id, gene_name, .before = dplyr::everything())

  # 保存全量 DEG
  all_csv <- file.path(out_dir_group, paste0("DEG_", group_col, ".csv"))
  readr::write_csv(markers, all_csv); message("Saved: ", all_csv)

  # 选每组 topN（按 avg_log2FC，先以显著性过滤）
  top_markers <- select_top_markers(markers, topn = TOPN_PER_GROUP, padj = PADJ_THRESH)

  # DotPlot（基于 topN；x 轴使用 gene_name）
  if (nrow(top_markers) > 0) {
    p_dot <- dotplot_from_markers(seu, top_markers, group_col, assay = assay)
    save_plot(p_dot, file.path(out_dir_group, paste0("DotPlot_", group_col, ".png")), DOT_WIDTH, DOT_HEIGHT)
    save_plot(p_dot, file.path(out_dir_group, paste0("DotPlot_", group_col, ".pdf")), DOT_WIDTH, DOT_HEIGHT)
  } else {
    message("No significant markers (adj.P < ", PADJ_THRESH, ") for ", group_col)
  }

  # 保存 topN 表
  top_csv <- file.path(out_dir_group, paste0("DEG_", group_col, "_top", TOPN_PER_GROUP, ".csv"))
  readr::write_csv(top_markers, top_csv); message("Saved: ", top_csv)

  # 每个 identity 的 topN FeaturePlot（标题=gene_name）
  feature_dir <- file.path(out_dir_group, "featureplots")
  featureplot_top_per_identity(seu, top_markers, group_col, feature_dir, assay = assay)
}

# ========= 主流程 =========
main <- function() {
  ensure_dir(OUTDIR)

  if (!file.exists(INPUT_RDS)) stop("INPUT_RDS not found: ", INPUT_RDS)
  message("Reading Seurat object: ", INPUT_RDS)
  seu <- readRDS(INPUT_RDS)

  if (!SOURCE_ASSAY %in% Assays(seu)) {
    stop("Assay '", SOURCE_ASSAY, "' not found in object. Available: ", paste(Assays(seu), collapse = ", "))
  }

  # 先基于 SOURCE_ASSAY 构建 RNA_sym（基于 gene_name 重命名行，复用 data 槽）
  seu <- build_RNA_sym(seu, source_assay = SOURCE_ASSAY, target_assay = ASSAY_USE,
                       gene_name_col = "gene_name", replace_underscore_with_dash = TRUE)

  # 元数据检查与组合列
  meta <- seu@meta.data
  need_cols <- c(AGE_COL, TYPE_COL)
  if (!all(need_cols %in% colnames(meta))) {
    stop("Missing metadata columns: ", paste(setdiff(need_cols, colnames(meta)), collapse = ", "))
  }
  meta[[COMBO_COL_NAME]] <- paste(meta[[AGE_COL]], meta[[TYPE_COL]], sep = "-")
  seu@meta.data <- meta

  # 准备 gene_id ↔ gene_name 映射（RNA_sym 已写入 meta.features）
  feature_map <- get_feature_map(seu, assay = ASSAY_USE)

  group_cols <- c(AGE_COL, TYPE_COL, COMBO_COL_NAME)

  for (gc in group_cols) {
    out_dir_group <- file.path(OUTDIR, safe_name(gc))
    ensure_dir(out_dir_group)

    # UMAP（与基因无关）
    umap_plot_one_group(seu, gc, out_dir_group)

    # DEG + DotPlot + FeaturePlot（全部在 RNA_sym 上，CSV 含 gene_id + gene_name）
    run_deg_one_group(seu, gc, out_dir_group, feature_map, assay = ASSAY_USE)
  }

  message("\nAll done. See outputs under: ", normalizePath(OUTDIR))
}

# ========= 命令行覆盖（可选）=========
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) INPUT_RDS <- args[1]
if (length(args) >= 2) OUTDIR    <- args[2]

main()

