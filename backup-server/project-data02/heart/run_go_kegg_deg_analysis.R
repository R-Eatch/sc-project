#!/usr/bin/env Rscript
# deg_enrich_go_kegg_from_deg_csv_BPCCMF_loop_clean.R  (rev-2025-07-30)
# -----------------------------------------------------------------------------
# 从按细胞类型与条件导出的 DEG CSV（示例：
#   DEG_Fibroblasts_up_in_DCM.csv, DEG_Fibroblasts_up_in_Donor.csv, ...）读取基因，
# 对每个 (celltype × condition) 的基因集合执行 GO/KEGG 富集：
#   * GO 采用 for 循环分别跑 BP / CC / MF，并各自保存 CSV 与绘图
#   * KEGG 可开关；自动 SYMBOL→ENTREZ 映射，并将 geneID 转回 SYMBOL 保存
# 输出一页组合 PDF：(BP | CC) / (MF | KEGG)
# 并把每个面板对应的表格写入 out_dir/tables/
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
})

# ===================== 全局配置 =====================
# 输入 DEG 目录（放置 csv）
deg_dir          <- "./deg_analysis"
# 输出目录（将创建 plots/ 与 tables/ 子目录）
out_dir          <- "./deg_enrich"
# 细胞类型与条件（用于文件名模板 `DEG_<celltype>_up_in_<condition>.csv`）
celltypes        <- c("Fibroblasts", "Cardiomyocytes", "Smooth_Muscle")
conditions       <- c("DCM", "Donor")
# 物种："human" 或 "mouse"
organism         <- "human"
# 是否执行 KEGG（TRUE/FALSE）
include_kegg     <- FALSE
# 基因筛选阈值（以 adj‑p 为主）
padj_filter      <- 0.01
logfc_filter     <- 0.1   # 仅保留 |logFC| >= 0.25 的基因；设为 0 可关闭该过滤
top_n_for_enrich <- 500    # 通过 |logFC| 排序后，选取用于 GO/KEGG 的前 N 个基因
min_gene_count   <- 3      # 用于富集的最小基因数
# 富集统计阈值
go_p_cutoff      <- 0.05
go_q_cutoff      <- 0.20
kegg_p_cutoff    <- 0.10
kegg_q_cutoff    <- 0.25
# 可视化参数
show_top         <- 15     # dotplot 展示的最多条目
plot_width       <- 14
plot_height      <- 12

# ===================== 辅助函数 =====================
messagef <- function(...) message(sprintf(...))

ensure_dirs <- function(base) {
  d1 <- file.path(base, "plots"); dir.create(d1, showWarnings = FALSE, recursive = TRUE)
  d2 <- file.path(base, "tables"); dir.create(d2, showWarnings = FALSE, recursive = TRUE)
  list(plots = d1, tables = d2)
}

ggsave_pdf_safe <- function(filename, plot, width, height) {
  ok <- TRUE
  tryCatch(
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf),
    error = function(e) {
      message("[warn] cairo_pdf not available, using default pdf: ", e$message)
      ok <<- FALSE
    }
  )
  if (!ok) {
    alt <- sub("\\.pdf$", "_fallback.pdf", filename, ignore.case = TRUE)
    tryCatch(
      ggsave(alt, plot, width = width, height = height),
      error = function(e2) message("[error] ggsave fallback failed: ", e2$message)
    )
  }
}

empty_plot <- function(txt) {
  ggplot() + annotate("text", x = .5, y = .5, label = txt, size = 4, lineheight = .9) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) + theme_void()
}

load_deg <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  names(df) <- tolower(names(df))
  if (!"gene" %in% names(df)) {
    if ("gene_name" %in% names(df)) df$gene <- df$gene_name else
    if ("genes"     %in% names(df)) df$gene <- df$genes else
      stop("DEG file lacks gene column (gene/gene_name/genes): ", path)
  }
  if (!"p_val_adj" %in% names(df) && !"p.val.adj" %in% names(df) && !"padj" %in% names(df)) {
    if ("p_val" %in% names(df)) df$p_val_adj <- p.adjust(df$p_val, method = "BH")
  }
  if (!"p_val_adj" %in% names(df)) {
    if ("p.val.adj" %in% names(df)) df$p_val_adj <- df$`p.val.adj` else
    if ("padj"       %in% names(df)) df$p_val_adj <- df$padj else
      stop("No adjusted-p column (p_val_adj/p.val.adj/padj), and p_val missing.")
  }
  df
}

species_setup <- function(which = c("mouse","human")) {
  which <- match.arg(which)
  if (which == "mouse") {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    list(OrgDb = org.Mm.eg.db, kegg = "mmu")
  } else {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    list(OrgDb = org.Hs.eg.db, kegg = "hsa")
  }
}

# 逐 ONTOLOGY 运行 GO，并返回绘图列表；同时写表
run_go_ontwise <- function(genes, OrgDb, table_dir, tag) {
  # 先准备占位图，后续仅在成功时替换
  plots <- list(
    BP = empty_plot("no GO BP"),
    CC = empty_plot("no GO CC"),
    MF = empty_plot("no GO MF")
  )
  for (ont in c("BP","CC","MF")) {
    go_res <- tryCatch(
      enrichGO(gene = genes, OrgDb = OrgDb, keyType = "SYMBOL", ont = ont,
               pAdjustMethod = "BH", pvalueCutoff = go_p_cutoff, qvalueCutoff = go_q_cutoff),
      error = function(e) { messagef("[warn] GO-%s failed for %s: %s", ont, tag, e$message); return(NULL) }
    )

    if (!is.null(go_res)) {
      res_df <- tryCatch(as.data.frame(go_res), error = function(e) NULL)
      if (!is.null(res_df) && nrow(res_df) > 0) {
        # 保存该本体的结果 CSV
        out_csv <- file.path(table_dir, sprintf("%s_GO_%s.csv", tag, ont))
        write.csv(res_df, out_csv, row.names = FALSE)

        # dotplot，showCategory 不超过可用条目数且至少为 1
        show_n <- max(1, min(show_top, nrow(res_df)))
        plots[[ont]] <- tryCatch(
          dotplot(go_res, showCategory = show_n) +
            ggtitle(sprintf("GO-%s", ont)) + theme_classic(base_size = 11),
          error = function(e) { messagef("[warn] GO-%s plot error for %s: %s", ont, tag, e$message);
                                empty_plot(sprintf("GO-%s plot error", ont)) }
        )
      }
    }
  }
  plots
}

run_kegg <- function(genes, OrgDb, kegg_org, table_dir, tag) {
  map <- suppressMessages(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb))
  if (is.null(map) || nrow(map) == 0) {
    messagef("[warn] %s : no ENTREZ mapping, skip KEGG", tag)
    return(list(plot = empty_plot("no KEGG (no map)"), res = NULL))
  }
  entrez <- unique(map$ENTREZID)
  kegg_res <- tryCatch(
    enrichKEGG(gene = entrez, organism = kegg_org, keyType = "ncbi-geneid",
               pAdjustMethod = "BH", pvalueCutoff = kegg_p_cutoff, qvalueCutoff = kegg_q_cutoff),
    error = function(e) { messagef("[warn] KEGG failed for %s: %s", tag, e$message); return(NULL) }
  )

  if (!is.null(kegg_res)) {
    kdf <- tryCatch(as.data.frame(kegg_res), error = function(e) NULL)
    if (!is.null(kdf) && nrow(kdf) > 0) {
      # Entrez → Symbol 便于阅读
      kdf$geneID <- vapply(kdf$geneID, function(x){
        ids <- strsplit(x, "/")[[1]]
        syms <- map$SYMBOL[match(ids, map$ENTREZID)]
        paste(na.omit(syms), collapse = "/")
      }, character(1))
      write.csv(kdf, file.path(table_dir, sprintf("%s_KEGG.csv", tag)), row.names = FALSE)

      show_n <- max(1, min(show_top, nrow(kdf)))
      kplot <- tryCatch(
        dotplot(kegg_res, showCategory = show_n) + ggtitle("KEGG") + theme_classic(base_size = 11),
        error = function(e) { messagef("[warn] KEGG plot error for %s: %s", tag, e$message);
                              empty_plot("KEGG plot error") }
      )
      return(list(plot = kplot, res = kegg_res))
    }
  }
  messagef("[info] %s : KEGG no significant terms", tag)
  list(plot = empty_plot("no KEGG"), res = NULL)
}

# ===================== 主流程 =====================
paths <- ensure_dirs(out_dir)
sp    <- species_setup(ifelse(tolower(organism) %in% c("human","homo","hs"), "human", "mouse"))
OrgDb <- sp$OrgDb; kegg_org <- sp$kegg

force_plot <- function(p, label="empty") {
  if (inherits(p, "ggplot") || inherits(p, "patchwork")) return(p)
  empty_plot(label)
}

for (ct in celltypes) {
  for (cond in conditions) {
    f <- file.path(deg_dir, sprintf("DEG_%s_up_in_%s.csv", ct, cond))
    if (!file.exists(f)) { messagef("[skip] missing: %s", f); next }

    df <- load_deg(f)

    # 检测 logFC 列（兼容 lfc / avg_log2FC / avg_logFC，已在 load_deg 中转为小写）
    lfc_col <- NULL
    if ("lfc"        %in% names(df)) lfc_col <- "lfc" else
    if ("avg_log2fc" %in% names(df)) lfc_col <- "avg_log2fc" else
    if ("avg_logfc"  %in% names(df)) lfc_col <- "avg_logfc"
    if (!is.null(lfc_col)) df$lfc_std <- df[[lfc_col]]

    # 先按 padj 过滤，再按 |logFC| 过滤
    df_padj <- dplyr::filter(df, p_val_adj < padj_filter)
    n_padj  <- nrow(df_padj)
    if (!is.null(lfc_col) && is.finite(logfc_filter) && logfc_filter > 0) {
      df_sel <- dplyr::filter(df_padj, abs(lfc_std) >= logfc_filter)
    } else if (is.null(lfc_col)) {
      messagef("[warn] %s | %s : logFC column not found; skip logFC filter", ct, cond)
      df_sel <- df_padj
    } else {
      # logfc_filter <= 0 表示不启用 logFC 过滤
      df_sel <- df_padj
    }

    # 基因排序并截取 TopN（若存在 logFC 列）
    if (!is.null(lfc_col)) {
      df_ranked <- df_sel %>% mutate(abs_lfc = abs(lfc_std)) %>% arrange(desc(abs_lfc))
      n_after <- nrow(df_ranked)
      genes <- unique(head(df_ranked$gene, top_n_for_enrich))
      messagef("[info] %s | %s : %d genes after padj<%.3g & |logFC|>=%.3g (from %d); taking top %d by |logFC| → %d genes",
               ct, cond, n_after, padj_filter, logfc_filter, n_padj, top_n_for_enrich, length(genes))
    } else {
      # 若无 logFC 列，无法按 |logFC| 排序与截断
      genes <- unique(df_sel$gene)
      messagef("[warn] %s | %s : logFC column not found; using all %d genes after padj/logFC filters (no topN)",
               ct, cond, length(genes))
    }

    if (length(genes) < min_gene_count) { messagef("[warn] %s | %s : < %d genes, skip", ct, cond, min_gene_count); next }

    tag <- sprintf("%s_%s", ct, cond)

    # GO：for 循环分别跑 BP/CC/MF，并保存各自 CSV
    go_plots <- run_go_ontwise(genes, OrgDb, paths$tables, tag)

    # KEGG：可选
    kegg_plot <- empty_plot("KEGG disabled")
    if (isTRUE(include_kegg)) {
      kp <- run_kegg(genes, OrgDb, kegg_org, paths$tables, tag)
      kegg_plot <- kp$plot
    }

    # 组合成一页 PDF（每格都是有效 ggplot）
    panel <- patchwork::wrap_plots(
      list(
        go_plots$BP,
        go_plots$CC,
        go_plots$MF,
        kegg_plot
      ), ncol = 2
    )
    pdf_path <- file.path(paths$plots, sprintf("%s_GO_KEGG.pdf", tag))
    ggsave_pdf_safe(pdf_path, panel, width = plot_width, height = plot_height)
    messagef("[ok] saved: %s", pdf_path)
  }
}

message("All enrichment tasks finished.")

