# ============================================================
# GO-BP enrichment per cell type (mouse) ─ RStudio version
# Requires: org.Mm.eg.db, clusterProfiler, dplyr, ggplot2, patchwork
# ============================================================

## 1. ── Parameters ------------------------------------------
input_csv   <- "conserved_gene_expression.csv"            # ← set your file path
plot_dir    <- "./result/plots"                   # PDFs go here
table_dir   <- "./result/go_tables"               # CSVs go here
show_n      <- 15                       # top terms to visualise

## 2. ── Libraries -------------------------------------------
suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(clusterProfiler)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

## 3. ── I/O prep --------------------------------------------
dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(c("celltype", "gene") %in% names(df))

## 4. ── Split by cell type ----------------------------------
gene_lists <- df %>%
  group_by(celltype) %>%
  summarise(genes = list(unique(gene)), .groups = "drop")

## 5. ── Loop over cell types --------------------------------
for (i in seq_len(nrow(gene_lists))) {
  ct    <- gene_lists$celltype[i]
  genes <- gene_lists$genes[[i]]
  
  message(sprintf("  %s  —  %d genes", ct, length(genes)))
  
  # 5b. GO-BP enrichment
  ego <- enrichGO(gene          = genes,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "SYMBOL",
                  ont           = "CC",
                  pAdjustMethod = "BH",
                  minGSSize     = 10,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  if (nrow(ego) == 0) {
    message(sprintf("  %s: no significant terms", ct))
    next
  }
  
  ## 5c. Save full result as CSV
  table_out <- file.path(table_dir, paste0(ct, "_GO_CC.csv"))
  write.csv(ego@result, table_out, row.names = FALSE)
  
  ## 5d. Prepare the top-15 subset for plotting
  ego_top <- ego %>% as.data.frame() %>% arrange(p.adjust) %>% head(show_n)
  # clusterProfiler’s plotting helpers need the enrichResult object,
  # so we rebuild a temporary object containing only the top terms:
  ego_top_obj <- ego
  ego_top_obj@result <- ego_top
  
  ## 5e. Plots
  p_bar <- barplot(ego_top_obj, showCategory = show_n) +
    ggtitle(ct) +
    theme_classic(base_size = 12)
  
  p_dot <- dotplot(ego_top_obj, showCategory = show_n) +
    ggtitle(ct) +
    theme_classic(base_size = 12)
  
  combined <- p_bar | p_dot
  
  ## 5f. Save PDF
  pdf_out <- file.path(plot_dir, paste0(ct, "_GO_CC.pdf"))
  ggsave(pdf_out, combined, width = 12, height = 6)
  message(sprintf(" — plots → %s, table → %s",
                  ct, basename(pdf_out), basename(table_out)))
}

message(" All done!")
