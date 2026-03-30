#!/usr/bin/env Rscript
# ==============================================================
# GO / KEGG enrichment per gene-set (dataset × type × celltype)
# - Input: 12 DEG CSVs (6 datasets × 2 types)
# - DEG CSV columns: group (celltype), names (gene SYMBOL), scores
# - For each celltype: take top 500 genes by scores and run enrichment
# - Output: tables/ and plots/ under result_root, filenames include dataset + type
# ==============================================================

## 1 ─ 用户设置 -------------------------------------------------
# 输入：DEG 结果根目录（按 type 分文件夹）
input_root   <- "D:/111/deg-result-imm"                # ← 含有 {type}/table/...
result_root  <- "D:/111/result-imm-enrich"             # ← 结果根目录

# 数据集与类型（共 12 份文件）
datasetlist <- list(
  'M-MG',
  'R-MG',
  'R-AG',
  'R-CG',
  'S-MG',
  'S-AG'
)
typelist <- c('celltype', 'newcelltype')

# 富集参数（保持你的流程逻辑不变）
ont          <- "MF"       # "BP" | "CC" | "MF" （仅用于 GO 分析）
output_go    <- TRUE
output_kegg  <- FALSE
show_n       <- 15
plot_format  <- "pdf"      # "pdf" or "png"

# DEG 取基因策略
TOP_N_GENES  <- 500           # ← 按 scores 排序后取前 500

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

## 4 ─ 主循环：逐文件（dataset × type）逐 celltype 分析 --------
for (type in typelist) {
  for (dataset in datasetlist) {
    
    input_csv <- file.path(input_root, type, "table",
                           sprintf("%s_%s_ranked_genes_for_annotation.csv", dataset, type))
    
    if (!file.exists(input_csv)) {
      message(sprintf("⚠️  Missing file, skipped: %s", input_csv))
      next
    }
    
    message("============================================================")
    message(sprintf("🧾  Reading: dataset=%s | type=%s", dataset, type))
    message(sprintf("     %s", input_csv))
    
    df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
    
    # 基本列检查（不改变流程，仅确保输入符合你描述的格式）
    required_cols <- c("group", "names", "scores")
    missing_cols  <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0) {
      message(sprintf("❌  File format error: missing columns [%s] in %s",
                      paste(missing_cols, collapse = ", "), input_csv))
      next
    }
    
    # 逐 celltype (group) 取 top500 基因列表
    celltypes <- unique(df$group)
    celltypes <- celltypes[!is.na(celltypes) & celltypes != ""]
    
    for (ct in celltypes) {
      
      sub_df <- df %>%
        filter(group == ct) %>%
        filter(!is.na(names), names != "") %>%
        arrange(desc(scores)) %>%
        slice_head(n = TOP_N_GENES)
      
      gene_vec <- sub_df$names %>% unique() %>% as.character()
      
      if (length(gene_vec) < 5) {
        message(sprintf("⚠️  %s | %s | %s skipped (<5 genes)", dataset, type, ct))
        next
      }
      
      message(sprintf("🟢  %s | %s | %s – %d genes (top %d by scores)",
                      dataset, type, ct, length(gene_vec), TOP_N_GENES))
      
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
                    file.path(table_dir, sprintf("%s_%s_%s_GO_%s.csv", dataset, type, ct, ont)),
                    row.names = FALSE)
          ## 绘图（top N）
          go_top        <- head(go_res[order(go_res$p.adjust), ], show_n)
          go_res@result <- go_top
          go_plot <- dotplot(go_res, showCategory = show_n) +
            ggtitle(sprintf("%s | %s | %s – GO-%s", dataset, type, ct, ont)) +
            theme_classic(base_size = 12)
          ggsave(file.path(plot_dir, sprintf("%s_%s_%s_GO_%s.%s", dataset, type, ct, ont, plot_format)),
                 go_plot, width = 7, height = 7)
        } else {
          message(sprintf("   ↳ %s | %s | %s: no significant GO terms", dataset, type, ct))
        }
      }
      
      # -------------------------------------------------- KEGG 分析（仍用 Entrez）
      if (output_kegg) {
        ent_map   <- bitr(gene_vec, fromType = "SYMBOL", toType = "ENTREZID",
                          OrgDb = org.Mm.eg.db)
        
        if (is.null(ent_map) || nrow(ent_map) == 0) {
          message(sprintf("   ↳ %s | %s | %s: no SYMBOL→ENTREZID mapping", dataset, type, ct))
          next
        }
        
        entrez_ids <- unique(ent_map$ENTREZID)
        if (length(entrez_ids) < 5) {
          message(sprintf("   ↳ %s | %s | %s: too few ENTREZ IDs after mapping (<5)", dataset, type, ct))
          next
        }
        
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
                    file.path(table_dir, sprintf("%s_%s_%s_KEGG.csv", dataset, type, ct)),
                    row.names = FALSE)
          ## 绘图
          kegg_top        <- head(kegg_res[order(kegg_res$p.adjust), ], show_n)
          kegg_res@result <- kegg_top
          kegg_plot <- dotplot(kegg_res, showCategory = show_n) +
            ggtitle(sprintf("%s | %s | %s – KEGG", dataset, type, ct)) +
            theme_classic(base_size = 12)
          ggsave(file.path(plot_dir, sprintf("%s_%s_%s_KEGG.%s", dataset, type, ct, plot_format)),
                 kegg_plot, width = 7, height = 7)
        } else {
          message(sprintf("   ↳ %s | %s | %s: no significant KEGG pathways", dataset, type, ct))
        }
      }
      
    } # end celltype loop
    
  } # end dataset loop
} # end type loop

message("🎉  Finished all dataset × type × celltype enrichments.")
