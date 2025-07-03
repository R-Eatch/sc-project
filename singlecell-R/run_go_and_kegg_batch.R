#!/usr/bin/env Rscript
# ==============================================================
# GO / KEGG enrichment per gene-set column   (mouse, SYMBOL for GO)
# ==============================================================
## 1 ─ 用户设置 -------------------------------------------------
input_csv    <- "D:/Rfile/result/0620/common_gene-0620.csv"                               # ← CSV 路径
result_root  <- "D:/Rfile/result/0620/"   # ← 结果根目录
ont          <- "MF"       # "BP" | "CC" | "MF"
output_go    <- TRUE
output_kegg  <- FALSE
show_n       <- 15
plot_format  <- "pdf"      # "pdf" or "png"

## 2 ─ 依赖包 ---------------------------------------------------
suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
})

## 3 ─ 输出目录 -------------------------------------------------
plot_dir  <- file.path(result_root, "plots")
table_dir <- file.path(result_root, "tables")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

## 4 ─ 读取数据 -------------------------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

## 5 ─ 主循环：逐列分析 ----------------------------------------
for (col in names(df)) {
  
  gene_vec <- df[[col]] %>% na.omit() %>% unique() %>% as.character()
  
  if (length(gene_vec) < 5) {
    message(sprintf("⚠️  %s skipped (<5 genes)", col))
    next
  }
  message(sprintf("🟢  %s – %d genes", col, length(gene_vec)))
  
  # -------------------------------------------------- GO 分析（SYMBOL）
  if (output_go) {
    go_res <- enrichGO(gene          = gene_vec,
                       OrgDb         = org.Mm.eg.db,
                       keyType       = "SYMBOL",   # ← 关键：直接用 SYMBOL
                       ont           = ont,
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.20)
    if (nrow(go_res)) {
      ## 保存表格
      write.csv(go_res@result,
                file.path(table_dir, sprintf("%s_GO_%s.csv", col, ont)),
                row.names = FALSE)
      ## 绘图（top N）
      go_top        <- head(go_res[order(go_res$p.adjust), ], show_n)
      go_res@result <- go_top
      go_plot <- dotplot(go_res, showCategory = show_n) +
        ggtitle(sprintf("%s – GO-%s", col, ont)) +
        theme_classic(base_size = 12)
      ggsave(file.path(plot_dir, sprintf("%s_GO_%s.%s", col, ont, plot_format)),
             go_plot, width = 7, height = 5)
    } else {
      message(sprintf("   ↳ %s: no significant GO terms", col))
    }
  }
  
  # -------------------------------------------------- KEGG 分析（仍用 Entrez）
  if (output_kegg) {
    ent_map   <- bitr(gene_vec, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Mm.eg.db)
    entrez_ids <- unique(ent_map$ENTREZID)
    
    kegg_res <- enrichKEGG(gene          = entrez_ids,
                           organism      = "mmu",
                           keyType       = "ncbi-geneid",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = 0.10)
    if (nrow(kegg_res)) {
      ## geneID 转回 SYMBOL
      kegg_res@result$geneID <- sapply(kegg_res@result$geneID, function(row) {
        ids  <- strsplit(row, "/")[[1]]
        syms <- ent_map$SYMBOL[match(ids, ent_map$ENTREZID)]
        paste(syms, collapse = "/")
      })
      ## 保存表格
      write.csv(kegg_res@result,
                file.path(table_dir, sprintf("%s_KEGG.csv", col)),
                row.names = FALSE)
      ## 绘图
      kegg_top        <- head(kegg_res[order(kegg_res$p.adjust), ], show_n)
      kegg_res@result <- kegg_top
      kegg_plot <- dotplot(kegg_res, showCategory = show_n) +
        ggtitle(sprintf("%s – KEGG", col)) +
        theme_classic(base_size = 12)
      ggsave(file.path(plot_dir, sprintf("%s_KEGG.%s", col, plot_format)),
             kegg_plot, width = 7, height = 5)
    } else {
      message(sprintf("   ↳ %s: no significant KEGG pathways", col))
    }
  }
}

message("🎉  Finished all columns.")
