library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(homologene)
library(patchwork)
library(tidyr)
#######################
# GO/KEGG enrichment pipeline (revised)
# ------------------------------------------------------------
# 1. The data‑reading, exporting, and enrichment calculation
#    logic is largely unchanged.
# 2. CSV files are first grouped by the **group** column; the
#    top *n* genes from each group are selected and merged into
#    one gene list that is used for enrichment.
# 3. One GO and one KEGG result are produced **per sample**
#    (e.g. M‑MG), not per group.
# 4. All exported files now include a global variable **lineage**
#    in their file names for easier bookkeeping.
# ------------------------------------------------------------

# ── libraries ───────────────────────────────────────────────
library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(homologene)
library(tidyr)

# ── global settings ─────────────────────────────────────────
lineage       <- "SEC"                    # tag added to every output file name
input_root    <- "D:/111/lineage-all-epi"     # where each <dataset>/<csv> lives
output_root   <- "D:/111/lineage-all-epi/lineage-GO-KEGG-result"          # where all results are written
ont           <- "BP"                         # ontology for GO enrichment
output_go     <- TRUE
output_kegg   <- TRUE
top_n         <- 30                           # top genes per group (n_gene)
# list of samples to iterate over
datasetlist <- c(#"M-MG", 
                 #"R-MG",
                 "S-MG"
                 #,"R-AG", "R-CG","S-AG"
                 )

# ── main loop ───────────────────────────────────────────────
for (dataset in datasetlist) {
  # 1. read csv ------------------------------------------------
  csv_path <- file.path(input_root, dataset, paste0(dataset, "_",lineage,"_pseudotime_genes.csv"))
  data     <- read.csv(csv_path)
  # 2. choose top_n genes per group, merge --------------------
  top_genes <- data %>%
    filter(logfoldchanges > 0) %>%
    group_by(group) %>%
    arrange(scores, .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup()
  
  genelist  <- unique(top_genes$names)
  
  # 3‑A. GO enrichment ----------------------------------------
  if (output_go) {
    ego <- enrichGO(
      gene          = genelist,
      OrgDb         = org.Mm.eg.db,
      keyType       = "SYMBOL",
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    go_plot <- dotplot(ego, showCategory = 15) +
      ggtitle(paste0(dataset, "_", lineage, "_GO"))
    
    # save plot
    ggsave(
      filename = paste0(dataset, "_", lineage, "_GO_Dotplot.png"),
      path     = output_root,
      plot     = go_plot,
      height   = 15, width = 10, limitsize = FALSE
    )
    ggsave(
      filename = paste0(dataset, "_", lineage, "_GO_Dotplot.pdf"),
      path     = output_root,
      plot     = go_plot,
      height   = 15, width = 10, limitsize = FALSE
    )
    
    # save table
    write.csv(
      ego@result,
      file.path(output_root, paste0("GO-result-", lineage, "-", dataset, ".csv")),
      row.names = FALSE
    )
    message(dataset, " GO enrichment finished")
  }
  
  # 3‑B. KEGG enrichment --------------------------------------
  if (output_kegg) {
    entrez_ids <- bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
    
    ekegg <- enrichKEGG(
      gene          = entrez_ids$ENTREZID,
      organism      = "mmu",
      keyType       = "ncbi-geneid",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.1
    )
    
    # replace ENTREZ IDs with gene symbols in geneID column
    if (nrow(ekegg@result) > 0) {
      ekegg@result$geneID <- sapply(ekegg@result$geneID, function(row) {
        ids  <- strsplit(row, "/")[[1]]
        syms <- entrez_ids$SYMBOL[match(ids, entrez_ids$ENTREZID)]
        paste(syms, collapse = "/")
      })
    }
    
    kegg_plot <- dotplot(ekegg, showCategory = 15) +
      ggtitle(paste0(dataset, "_", lineage, "_KEGG"))
    
    # save plot
    ggsave(
      filename = paste0(dataset, "_", lineage, "_KEGG_Dotplot.png"),
      path     = output_root,
      plot     = kegg_plot,
      height   = 15, width = 10, limitsize = FALSE
    )
    ggsave(
      filename = paste0(dataset, "_", lineage, "_KEGG_Dotplot.pdf"),
      path     = output_root,
      plot     = kegg_plot,
      height   = 15, width = 10, limitsize = FALSE
    )
    
    # save table
    write.csv(
      ekegg@result,
      file.path(output_root, paste0("KEGG-result-", lineage, "-", dataset, ".csv")),
      row.names = FALSE
    )
    message(dataset, " KEGG enrichment finished")
  }
}
