#!/usr/bin/env Rscript
# ==============================================================
# One-column Top-N enrichment (GO-BP/CC/MF + KEGG) → single PDF
# ==============================================================

## 0 ─ 用户参数 --------------------------------------------------
input_csv    <- "D:/111/LumHR_unique_results/LummHR-SCGB_human_unique_inData.csv"
gene_col     <- "mouse_gene"    # 基因名列 (SYMBOL)
rank_col     <- "scores"        # 排序列
top_n_genes  <- 300             # Top-N
dataset      <- "LummHR-SCGB"        # ←★ 新增：数据集标签 (影响文件名 & 图标题)
result_root  <- "D:/111/LumHR_unique_results/enrich_result"

plot_width   <- 20              # ←★ 放大宽度
plot_height  <- 20              # ←★ 放大高度
## --------------------------------------------------------------

## 1 ─ 依赖包 ----------------------------------------------------
suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## 2 ─ 输出目录 --------------------------------------------------
plot_dir  <- file.path(result_root, "plots")
table_dir <- file.path(result_root, "tables")
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

## 3 ─ 读取并选 Top-N 基因 --------------------------------------
df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(c(gene_col, rank_col) %in% names(df)))
  stop("gene_col 或 rank_col 不存在于数据框。")

top_genes <- df %>% 
  arrange(desc(.data[[rank_col]])) %>% 
  pull(.data[[gene_col]]) %>% 
  na.omit() %>% unique() %>% head(top_n_genes)

cat(sprintf("🟢  取到 %d 个基因用于富集分析\n", length(top_genes)))

## 4 ─ 占位图 & 写表辅助 ---------------------------------------
empty_plot <- function(title_txt) {
  ggplot() +
    annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
             fill = NA, colour = "grey30", linetype = "dashed") +
    annotate("text", x = .5, y = .5, label = title_txt, size = 4) +
    theme_void()
}
save_table <- function(res_obj, path) {
  write.csv(as.data.frame(res_obj), path, row.names = FALSE)
}

## 5 ─ GO (BP/CC/MF) --------------------------------------------
go_onts  <- c("BP","CC","MF")
go_plots <- list()

for (ont in go_onts) {
  go_res <- enrichGO(
    gene          = top_genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.20
  )
  if (nrow(go_res)) {
    save_table(
      go_res@result,
      file.path(table_dir, sprintf("%s_GO_%s.csv", dataset, ont))
    )
    go_top        <- head(go_res[order(go_res$p.adjust), ], 15)
    go_res@result <- go_top
    go_plots[[ont]] <- dotplot(go_res, showCategory = 15) +
      ggtitle(sprintf("%s – GO-%s", dataset, ont)) +
      theme_classic(base_size = 11)
  } else {
    go_plots[[ont]] <- empty_plot(sprintf("%s: no GO-%s terms", dataset, ont))
  }
}

## 6 ─ KEGG ------------------------------------------------------
entrez_map <- bitr(top_genes, fromType = "SYMBOL", toType = "ENTREZID",
                   OrgDb = org.Mm.eg.db)
entrez_ids <- unique(entrez_map$ENTREZID)

kegg_res <- enrichKEGG(gene = entrez_ids,
                       organism = "mmu",
                       keyType = "ncbi-geneid",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.10)

if (nrow(kegg_res)) {
  kegg_res@result$geneID <- sapply(kegg_res@result$geneID, function(row) {
    ids  <- strsplit(row, "/")[[1]]
    syms <- entrez_map$SYMBOL[match(ids, entrez_map$ENTREZID)]
    paste(syms, collapse = "/")
  })
  save_table(kegg_res@result,
             file.path(table_dir, sprintf("%s_KEGG.csv", dataset)))
  
  kegg_top        <- head(kegg_res[order(kegg_res$p.adjust), ], 15)
  kegg_res@result <- kegg_top
  
  kegg_plot <- dotplot(kegg_res, showCategory = 15) +
    ggtitle(sprintf("%s – KEGG", dataset)) +
    theme_classic(base_size = 11)
} else {
  kegg_plot <- empty_plot(sprintf("%s: no KEGG pathways", dataset))
}

## 7 ─ 拼图并输出 PDF -------------------------------------------
combined <- (go_plots[["BP"]] | go_plots[["CC"]]) /
  (go_plots[["MF"]] | kegg_plot)

pdf(file.path(plot_dir,
              sprintf("Combined_GO_KEKG_%s.pdf", dataset)),
    width = plot_width, height = plot_height)
print(combined)
dev.off()

cat("🎉  All done! 结果与图表已输出至：", result_root, "\n")
