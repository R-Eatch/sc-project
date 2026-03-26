#!/usr/bin/env Rscript
# quick_subset_and_deg.R   (FINAL rev-2025-07-25, data-mask fix)
# -----------------------------------------------------------------------------
# 功能：
#   1. 读取 Seurat v4 .Robj.gz → 用 RNA counts 重建 v5 对象。
#   2. 对 *完整* 数据集执行标准流程：LogNormalize → HVG → Scale → PCA → Harmony → UMAP。
#   3. 绘制并保存描绘 *完整* 数据集的 4-in-1 UMAP 对比图。
#   4. 按 meta.data$Names 筛选保留 Fibroblasts / Cardiomyocytes / Smooth_Muscle。
#   5. 对筛选出的 *每种* 细胞类型，按 Condition 进行差异表达基因 (DEG) 分析。
#   6. 为每种细胞类型的 DEG 结果保存 CSV 表格和 DotPlot 可视化图。
#   7. 将最终的子集对象按 Condition 拆分并保存为 *.rds 文件，并额外保存未分 condition 的整体子集对象。
# -----------------------------------------------------------------------------
# 用法：
#   Rscript quick_subset_and_deg.R
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)          # v5
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(harmony)         # >=1.0.0
  library(patchwork)       # for panel layout
})
options(future.globals.maxSize = 10 * 1024^3)
## ---------------------------- 配置 -----------------------------------------
# --- 输入/输出文件 ---
file_in         <- "./data/GSE183852_DCM_Nuclei.Robj.gz"
output_prefix   <- "./data/heart_Fib_Cardio_Smooth"      # <prefix>_<Condition>.rds
plot_out_pdf    <- "umap_compare_FULL.pdf"
plot_out_png    <- "umap_compare_FULL.png"
deg_out_dir     <- "./deg_analysis"                      # DEG 结果输出目录

# --- 细胞类型筛选 ---
celltype_col    <- "Names"
celltypes_keep  <- c("Fibroblasts", "Cardiomyocytes", "Smooth_Muscle")

# --- PCA / Harmony / UMAP 设置 ---
n_pcs_init      <- 50   # 先计算 50 PCs
n_pcs_use       <- 10   # 用前 10 PCs 做 UMAP / neighbor
k_neighbors     <- 20
resolution_use  <- 0.2
set.seed(42)

# --- DEG 分析设置 ---
deg_logfc_threshold <- 0.25
deg_min_pct         <- 0.1
deg_top_n_dotplot   <- 20   # 在 DotPlot 中显示 Top N 基因（本版每个 condition Top10 合并）

## ---------------------------------------------------------------------------

# 创建输出目录
if (!dir.exists(deg_out_dir)) {
  dir.create(deg_out_dir, recursive = TRUE)
}

logmsg <- function(...) cat("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n", sep = "")

## 1. 读取并重建对象 ---------------------------------------------------------
logmsg("Loading:", file_in)
load_env <- new.env()
obj_name <- load(file_in, envir = load_env)[1]
old_obj  <- get(obj_name, envir = load_env)
logmsg("Loaded v4 object:", obj_name)

# 仅用 RNA counts
counts <- if ("RNA" %in% names(old_obj@assays)) {
  logmsg("Using counts from 'RNA' assay")
  old_obj@assays$RNA@counts
} else {
  first_assay <- names(old_obj@assays)[1]
  logmsg("'RNA' assay not found, using counts from", first_assay)
  old_obj@assays[[first_assay]]@counts
}
meta <- old_obj@meta.data

seu <- CreateSeuratObject(counts = counts, meta.data = meta, project = "SeuratProject")
logmsg("v5 object created | Cells:", ncol(seu), " Genes:", nrow(seu))

# 复制原 UMAP 如存在
if ("umap" %in% Reductions(old_obj)) {
  logmsg("Copying original UMAP to 'umap_orig'")
  emb <- Embeddings(old_obj, "umap")
  seu[["umap_orig"]] <- CreateDimReducObject(embeddings = emb, key = "origUMAP_")
} else {
  # 如果不存在原始UMAP，创建一个空的reduction以便后续绘图代码统一
  logmsg("Original UMAP not found. Skipping related plot.")
  p1 <- ggplot() + theme_void() + ggtitle("Orig UMAP • Names (Not Found)")
  p2 <- ggplot() + theme_void() + ggtitle("Orig UMAP • Condition (Not Found)")
}

## 2. 对 *完整* 数据集进行预处理 ---------------------------------------------
logmsg("--- Pre-processing full dataset ---")

logmsg("LogNormalize (scale.factor=1e4)…")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)

logmsg("Finding 3k HVGs (vst)…")
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

logmsg("Scaling all genes…")
seu <- ScaleData(seu, features = rownames(seu), verbose = FALSE)

logmsg("Running PCA (", n_pcs_init, " PCs)…")
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = n_pcs_init, verbose = FALSE)

logmsg("Running Harmony on 'orig.ident'…")
seu <- RunHarmony(
  object        = seu,
  group.by.vars = "orig.ident"
)

logmsg("Running UMAP (dims 1-", n_pcs_use, ")…")
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:n_pcs_use, verbose = FALSE)

logmsg("Finding neighbors & clusters…")
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:n_pcs_use, k.param = k_neighbors, verbose = FALSE)
seu <- FindClusters(seu, resolution = resolution_use, verbose = FALSE)

## 3. 绘制并保存 *完整* 数据集的 UMAP 对比图 ---------------------------------
logmsg("Generating comparison plots for the full dataset…")
if (exists("p1")) { # 如果原始UMAP不存在，则使用占位符图
  # p1, p2 already created
} else {
  p1 <- DimPlot(seu, reduction = "umap_orig", group.by = "Names", pt.size = 0.3) + ggtitle("Orig UMAP • Names")
  p2 <- DimPlot(seu, reduction = "umap_orig", group.by = "Condition", pt.size = 0.3) + ggtitle("Orig UMAP • Condition")
}
p3 <- DimPlot(seu, reduction = "umap", group.by = "Names", pt.size = 0.3) + ggtitle("New UMAP • Names")
p4 <- DimPlot(seu, reduction = "umap", group.by = "Condition", pt.size = 0.3) + ggtitle("New UMAP • Condition")
combo <- (p1 | p2) / (p3 | p4) + plot_layout(guides = "collect")

ggsave(plot_out_pdf, combo, width = 12, height = 10)
ggsave(plot_out_png, combo, width = 12, height = 10, dpi = 300)
logmsg("Saved", plot_out_pdf, "and", plot_out_png)

## 4. 子集保留指定细胞类型 ----------------------------------------------------
if (!celltype_col %in% colnames(seu@meta.data))
  stop("Column '", celltype_col, "' not found in meta.data")

logmsg("--- Subsetting for target cell types ---")
logmsg("Subsetting cell types:", paste(celltypes_keep, collapse = ", "))
seu_sub <- subset(seu, subset = (!!sym(celltype_col)) %in% celltypes_keep)
# 确保因子水平被移除，避免后续步骤出错
seu_sub@meta.data[[celltype_col]] <- droplevels(seu_sub@meta.data[[celltype_col]])
logmsg("Cells after filtering:", ncol(seu_sub))

## 5. 对每种细胞类型进行 DEG 分析并可视化 -------------------------------------
logmsg("--- Running DEG analysis for each cell type ---")
if (!"Condition" %in% colnames(seu_sub@meta.data))
  stop("meta.data does not contain 'Condition' column")

conditions <- unique(seu_sub@meta.data$Condition)
if (length(conditions) != 2) {
  logmsg("Warning: Expected 2 conditions for DEG analysis, but found", length(conditions),
         ". Skipping DEG analysis. Conditions found:", paste(conditions, collapse=", "))
} else {
  # 固定比较顺序以确保方向一致：正 logFC = 在 DCM 上调，负 logFC = 在 Donor 上调
  target_order <- c("Donor", "DCM")
  for (ct in celltypes_keep) {
    logmsg("Processing:", ct)
    ct_sub <- subset(seu_sub, subset = (!!sym(celltype_col)) == ct)

    # 该细胞类型是否在两个条件下都存在？
    present_conds <- unique(ct_sub$Condition)
    if (length(intersect(target_order, present_conds)) < 2) {
      logmsg("  -> Skipping DEG for", ct, "- not present in both conditions.")
      next
    }

    # 设定恒定的因子顺序与比较方向
    ct_sub$Condition <- factor(ct_sub$Condition, levels = target_order)
    ident_1 <- "DCM"   # 病例组
    ident_2 <- "Donor" # 对照组

    Idents(ct_sub) <- "Condition"
    logmsg("  -> Running FindMarkers for", ct, "between", ident_1, "and", ident_2)

    deg_results <- FindMarkers(
      ct_sub,
      ident.1         = ident_1,
      ident.2         = ident_2,
      logfc.threshold = deg_logfc_threshold,
      min.pct         = deg_min_pct,
      verbose         = FALSE
    )

    # 识别 logFC 列名（兼容 Seurat 版本），并创建显式列 `lfc` 以避免 data mask 问题
    lfc_col <- if ("avg_log2FC" %in% colnames(deg_results)) "avg_log2FC" else
               if ("avg_logFC"   %in% colnames(deg_results))   "avg_logFC"   else
               stop("DEG 结果中未找到 avg_log2FC/avg_logFC 列。")

    # 若缺少 p_val_adj，则由 p_val 生成，以供排序
    if (!"p_val_adj" %in% colnames(deg_results) && "p_val" %in% colnames(deg_results)) {
      deg_results$p_val_adj <- p.adjust(deg_results$p_val, method = "BH")
    }

    # 注释方向、复制列，避免 .data 取列导致的 data mask 报错
    deg_results$gene <- rownames(deg_results)
    deg_results$lfc  <- deg_results[[lfc_col]]
    deg_results$abs_logFC <- abs(deg_results$lfc)
    deg_results$up_in     <- ifelse(deg_results$lfc > 0, ident_1, ident_2)
    deg_results$direction <- ifelse(deg_results$lfc > 0, paste0("Up in ", ident_1),
                                    ifelse(deg_results$lfc < 0, paste0("Up in ", ident_2), "No change"))
    deg_results <- deg_results[order(deg_results$p_val_adj, -deg_results$abs_logFC), , drop = FALSE]

    # 保存总表与方向分表
    deg_csv_file <- file.path(deg_out_dir, sprintf("DEG_%s_%s_vs_%s.csv", ct, ident_1, ident_2))
    write.csv(deg_results, deg_csv_file, row.names = FALSE)
    logmsg("  -> Saved DEG results to ", deg_csv_file)

    deg_csv_up_dcm   <- file.path(deg_out_dir, sprintf("DEG_%s_up_in_%s.csv", ct, ident_1))
    deg_csv_up_donor <- file.path(deg_out_dir, sprintf("DEG_%s_up_in_%s.csv", ct, ident_2))
    write.csv(deg_results[deg_results$lfc > 0, , drop = FALSE], deg_csv_up_dcm,   row.names = FALSE)
    write.csv(deg_results[deg_results$lfc < 0, , drop = FALSE], deg_csv_up_donor, row.names = FALSE)

    # DotPlot：每个 condition 各取 Top10（先按 adj‑p，再按 |logFC|），合并绘制一张
    if (nrow(deg_results) > 0) {
      n_per_cond <- 10L
      pos_genes <- head(deg_results[deg_results$lfc > 0, "gene", drop = TRUE], n_per_cond)
      neg_genes <- head(deg_results[deg_results$lfc < 0, "gene", drop = TRUE], n_per_cond)
      top_genes <- unique(c(pos_genes, neg_genes))

      if (length(top_genes) > 0) {
        subtitle_txt <- paste0("Top10 per condition • ", ident_1, "=", length(pos_genes),
                               ", ", ident_2, "=", length(neg_genes))
        dotplot <- DotPlot(ct_sub, features = top_genes, group.by = "Condition") +
          ggtitle(paste("DEGs in", ct), subtitle = subtitle_txt) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))

        plot_pdf <- file.path(deg_out_dir, sprintf("DotPlot_DEG_%s.pdf", ct))
        plot_png <- file.path(deg_out_dir, sprintf("DotPlot_DEG_%s.png", ct))

        ggsave(plot_pdf, dotplot, width = 10, height = 6)
        ggsave(plot_png, dotplot, width = 10, height = 6, dpi = 300)
        logmsg("  -> Saved DotPlot to ", plot_pdf, " and ", plot_png)
      } else {
        logmsg("  -> No genes passed thresholds for DotPlot in ", ct)
      }
    } else {
      logmsg("  -> No significant DEGs found for ", ct, " with current thresholds. Skipping DotPlot.")
    }
  }
}

## 6. 拆分并保存最终的 RDS 对象 -----------------------------------------------
logmsg("--- Splitting and saving final subset objects ---")
final_conds <- unique(seu_sub@meta.data$Condition)
logmsg("Conditions in the final subset:", paste(final_conds, collapse = ", "))

# 额外：保存未分 condition 的整体子集对象
fout_all <- sprintf("%s_ALL.rds", output_prefix)
logmsg("Saving ", fout_all)
saveRDS(seu_sub, fout_all)

for (cond in final_conds) {
  sub_to_save <- subset(seu_sub, subset = Condition == cond)
  fout <- sprintf("%s_%s.rds", output_prefix, cond)
  logmsg("Saving ", fout)
  saveRDS(sub_to_save, fout)
}

logmsg("All done ✓")

