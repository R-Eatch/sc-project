#!/usr/bin/env Rscript
# deg_enrich_by_cluster.R  (rev-2025-09-29)
# 读取 ./result/{group_name}/DEG_{group_name}.csv
# 按 CSV 中的 `cluster` 列分组；以 `gene_name`(或等价列) 作为基因符号做 GO/KEGG。
# 为每个 (group_name × cluster) 产出：
#   - tables/: GO-BP/CC/MF 与 KEGG 的结果 CSV
#   - plots/: 1 页组合 PDF：BP | CC
#                                 MF | KEGG
# 物种固定 mouse。

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)  # mouse
})

# ===================== 配置 =====================
# 输入根目录与分组名称
result_root       <- "./result"
group_names       <- c("age_celltype", "Age_group", "Main_cell_type")

# 输出根目录（每个 group 会在其各自目录下新建 enrich_* 子目录）
enrich_dir_name   <- "enrich_by_cluster"

# 富集与筛选参数（与你原脚本一致/相近，可按需修改）
include_kegg      <- TRUE         # 是否做 KEGG
padj_filter       <- 0.01         # 仅保留显著基因
logfc_filter      <- 0.10         # 仅保留 |logFC| >= 0.10（<=0 则不启用该过滤）
top_n_for_enrich  <- 500          # 取按 |logFC| 排序的 Top N
min_gene_count    <- 3            # 富集最少基因数

go_p_cutoff       <- 0.05
go_q_cutoff       <- 0.20
kegg_p_cutoff     <- 0.10
kegg_q_cutoff     <- 0.25

show_top          <- 15           # dotplot 展示条目上限
plot_width        <- 14
plot_height       <- 12

# ===================== 工具函数 =====================
messagef <- function(...) message(sprintf(...))

ensure_dirs <- function(base) {
  d1 <- file.path(base, "plots");  dir.create(d1, showWarnings = FALSE, recursive = TRUE)
  d2 <- file.path(base, "tables"); dir.create(d2, showWarnings = FALSE, recursive = TRUE)
  list(plots = d1, tables = d2)
}

ggsave_pdf_safe <- function(filename, plot, width, height) {
  ok <- TRUE
  tryCatch(
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf),
    error = function(e) { message("[warn] cairo_pdf not available: ", e$message); ok <<- FALSE }
  )
  if (!ok) {
    alt <- sub("\\.pdf$", "_fallback.pdf", filename, ignore.case = TRUE)
    tryCatch(ggsave(alt, plot, width = width, height = height),
             error = function(e2) message("[error] ggsave fallback failed: ", e2$message))
  }
}

empty_plot <- function(txt) {
  ggplot() + annotate("text", x=.5, y=.5, label=txt, size=4, lineheight=.9) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) + theme_void()
}

# 统一读取并标准化列名；确保 gene / p_val_adj / cluster 等存在
load_deg <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  names(df) <- tolower(names(df))

  # gene 列：优先 gene_name，其次 gene / genes / gene_nam（防 OCR/截断）
  if (!"gene_name" %in% names(df)) {
    if ("gene"     %in% names(df)) df$gene_name <- df$gene else
    if ("genes"    %in% names(df)) df$gene_name <- df$genes else
    if ("gene_nam" %in% names(df)) df$gene_name <- df$gene_nam
  }
  if (!"gene_name" %in% names(df)) stop("No gene_name/gene/genes/gene_nam column in: ", path)

  # p_val_adj：若无则由 p_val 计算
  if (!"p_val_adj" %in% names(df)) {
    if ("padj" %in% names(df)) df$p_val_adj <- df$padj else
    if ("p.val.adj" %in% names(df)) df$p_val_adj <- df$`p.val.adj` else
    if ("p_val" %in% names(df)) df$p_val_adj <- p.adjust(df$p_val, method = "BH")
  }
  if (!"p_val_adj" %in% names(df)) stop("No adjusted-p column (p_val_adj/padj/p.val.adj) and no p_val to compute.")

  # logFC 统一到 lfc_std，兼容 avg_log2fc / avg_logfc / lfc
  lfc_col <- NULL
  if ("avg_log2fc" %in% names(df)) lfc_col <- "avg_log2fc" else
  if ("avg_logfc"  %in% names(df)) lfc_col <- "avg_logfc"  else
  if ("lfc"        %in% names(df)) lfc_col <- "lfc"
  if (!is.null(lfc_col)) df$lfc_std <- df[[lfc_col]]

  # cluster 列
  if (!"cluster" %in% names(df)) stop("No 'cluster' column in: ", path)

  df
}

run_go_ontwise <- function(genes, OrgDb, table_dir, tag) {
  plots <- list(BP=empty_plot("no GO BP"), CC=empty_plot("no GO CC"), MF=empty_plot("no GO MF"))
  for (ont in c("BP","CC","MF")) {
    go_res <- tryCatch(
      enrichGO(gene=genes, OrgDb=OrgDb, keyType="SYMBOL", ont=ont,
               pAdjustMethod="BH", pvalueCutoff=go_p_cutoff, qvalueCutoff=go_q_cutoff),
      error=function(e){ messagef("[warn] GO-%s failed for %s: %s", ont, tag, e$message); NULL }
    )
    if (!is.null(go_res)) {
      res_df <- tryCatch(as.data.frame(go_res), error=function(e) NULL)
      if (!is.null(res_df) && nrow(res_df) > 0) {
        write.csv(res_df, file.path(table_dir, sprintf("%s_GO_%s.csv", tag, ont)), row.names=FALSE)
        show_n <- max(1, min(show_top, nrow(res_df)))
        plots[[ont]] <- tryCatch(
          dotplot(go_res, showCategory=show_n) + ggtitle(sprintf("GO-%s", ont)) + theme_classic(base_size = 11),
          error=function(e){ messagef("[warn] GO-%s plot error for %s: %s", ont, tag, e$message);
                             empty_plot(sprintf("GO-%s plot error", ont)) }
        )
      }
    }
  }
  plots
}

run_kegg <- function(genes, OrgDb, kegg_org, table_dir, tag) {
  map <- suppressMessages(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb))
  if (is.null(map) || nrow(map) == 0) return(list(plot=empty_plot("no KEGG (no map)"), res=NULL))
  entrez <- unique(map$ENTREZID)
  kegg_res <- tryCatch(
    enrichKEGG(gene=entrez, organism=kegg_org, keyType="ncbi-geneid",
               pAdjustMethod="BH", pvalueCutoff=kegg_p_cutoff, qvalueCutoff=kegg_q_cutoff),
    error=function(e){ messagef("[warn] KEGG failed for %s: %s", tag, e$message); NULL }
  )
  if (!is.null(kegg_res)) {
    kdf <- tryCatch(as.data.frame(kegg_res), error=function(e) NULL)
    if (!is.null(kdf) && nrow(kdf) > 0) {
      kdf$geneID <- vapply(kdf$geneID, function(x){
        ids <- strsplit(x, "/")[[1]]
        syms <- map$SYMBOL[match(ids, map$ENTREZID)]
        paste(na.omit(syms), collapse = "/")
      }, character(1))
      write.csv(kdf, file.path(table_dir, sprintf("%s_KEGG.csv", tag)), row.names = FALSE)
      show_n <- max(1, min(show_top, nrow(kdf)))
      kplot <- tryCatch(
        dotplot(kegg_res, showCategory=show_n) + ggtitle("KEGG") + theme_classic(base_size = 11),
        error=function(e){ messagef("[warn] KEGG plot error for %s: %s", tag, e$message); empty_plot("KEGG plot error") }
      )
      return(list(plot=kplot, res=kegg_res))
    }
  }
  list(plot=empty_plot("no KEGG"), res=NULL)
}

# ===================== 主流程 =====================
OrgDb    <- org.Mm.eg.db
kegg_org <- "mmu"

for (grp in group_names) {
  deg_csv  <- file.path(result_root, grp, sprintf("DEG_%s.csv", grp))
  if (!file.exists(deg_csv)) { messagef("[skip] %s not found", deg_csv); next }

  df <- load_deg(deg_csv)

  # 输出目录：./result/<group>/enrich_by_cluster/{plots,tables}
  out_base <- file.path(result_root, grp, enrich_dir_name)
  paths    <- ensure_dirs(out_base)

  # 按 cluster 分组
  clusters <- sort(unique(df$cluster))
  messagef("[info] Group %s: %d clusters detected", grp, length(clusters))

  for (cl in clusters) {
    sub <- df %>% filter(cluster == cl)

    # padj 过滤
    sub <- sub %>% filter(p_val_adj < padj_filter)
    if (nrow(sub) == 0) { messagef("[warn] %s | cluster=%s : no genes after padj<%.3g", grp, cl, padj_filter); next }

    # logFC 过滤（如存在）
    if ("lfc_std" %in% names(sub) && is.finite(logfc_filter) && logfc_filter > 0) {
      sub <- sub %>% filter(abs(lfc_std) >= logfc_filter)
      if (nrow(sub) == 0) { messagef("[warn] %s | %s : no genes after |logFC|>=%.3g", grp, cl, logfc_filter); next }
    } else if (!("lfc_std" %in% names(sub))) {
      messagef("[warn] %s | %s : logFC column not found; skip logFC filter", grp, cl)
    }

    # 排序并取 TopN（若有 logFC）
    genes <- NULL
    if ("lfc_std" %in% names(sub)) {
      sub <- sub %>% mutate(abs_lfc = abs(lfc_std)) %>% arrange(desc(abs_lfc))
      genes <- unique(head(sub$gene_name, top_n_for_enrich))
    } else {
      genes <- unique(sub$gene_name)
    }

    # 去掉 NA/空
    genes <- unique(na.omit(trimws(genes)))
    if (length(genes) < min_gene_count) {
      messagef("[warn] %s | %s : < %d genes, skip", grp, cl, min_gene_count); next
    }

    tag <- sprintf("%s__%s", grp, gsub("[/\\s]+","_", cl))

    # GO
    go_plots <- run_go_ontwise(genes, OrgDb, paths$tables, tag)

    # KEGG
    kegg_plot <- empty_plot("KEGG disabled")
    if (isTRUE(include_kegg)) {
      kp <- run_kegg(genes, OrgDb, kegg_org, paths$tables, tag)
      kegg_plot <- kp$plot
    }

    # 组合面板并保存
    panel <- patchwork::wrap_plots(
      list(go_plots$BP, go_plots$CC, go_plots$MF, kegg_plot),
      ncol = 2
    )
    pdf_path <- file.path(paths$plots, sprintf("%s_GO_KEGG.pdf", tag))
    ggsave_pdf_safe(pdf_path, panel, width = plot_width, height = plot_height)
    messagef("[ok] saved: %s", pdf_path)
  }
}

message("All enrichment tasks finished.")

