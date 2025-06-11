# Monocle3 pseudotime analysis with automated root selection via principal graph nodes

library(Seurat)
library(monocle3)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)

#### 定义全局变量 ###
dataset <- "S-MAG"
cores <- 1
####################

##### 读取数据并构建 CDS #####
ea <- readRDS(paste0(dataset, "_for_monocle.rds"))
data <- GetAssayData(ea, assay = 'RNA', slot = 'counts')
cell_metadata <- ea@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
rm(data, cell_metadata, gene_annotation)

# 预处理与降维
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, cores = cores)
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "newcelltype")

# 使用 Seurat UMAP embedding
int.embed <- Embeddings(ea, reduction = "umap")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "newcelltype")
save_monocle_objects(cds=cds, directory_path=paste0(dataset,'_cds_for_subset.rds'))
#############电脑上操作#############
cds_subset <- load_monocle_objects(directory_path=paste0("D:/111/",dataset,'_cds_for_subset.rds'))


#cds_subset <- choose_cells(cds)  ####手动操作####
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset,use_partition = FALSE)
cds_subset <- order_cells(cds_subset)   ####手动操作####

save_monocle_objects(cds=cds_subset, directory_path=paste0(dataset,'_subset_cds.rds'))
#######################################
cds[]
#########绘图##################

cds_subset <- load_monocle_objects(directory_path=paste0(dataset,'_subset_cds.rds'))
p1 <- plot_cells(cds_subset,
                 color_cells_by = "pseudotime",
                 label_cell_groups=TRUE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 graph_label_size=3)
print(p1)
p2 <- plot_cells(cds_subset,
                 color_cells_by = "newcelltype",
                 label_groups_by_cluster=TRUE,
                 label_leaves=TRUE,
                 label_branch_points=TRUE,
                 graph_label_size=3)
print(p2)

# 使用plot_grid来组合图表，设置一行三个图
p1p2 <- wrap_plots(p1,p2)

# 保存featureplot，调整width和height以适应新的布局
ggsave(paste0(dataset, "_MONOCLE3.png"), p1p2, width = 16, height = 8)
#######################
#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: find_branch_genes_S-MAG.R
# Purpose: Identify branch‑dependent genes distinguishing Mammary Gland (MG) vs
#          Apocrine Gland (AG) developmental fates in the sugar‑glider epithelium
# Author : Adrian  (generated via ChatGPT‑o3)
# ------------------------------------------------------------------------------
# Expected input
#   •   S‑MAG_cds_rooted.rds  – a Monocle3 cell_data_set that already has:
#           - raw counts matrix (assay_table)
#           - UMAP reduction + principal graph (learn_graph)
#           - pseudotime values (order_cells)
#           - colData columns:  newcelltype, stage, gland ("MG"/"AG")
# ------------------------------------------------------------------------------
# Output files (written to working directory)
#   1.  MG_vs_AG_branch_genes_q0.05.csv          – ranked list of branch genes
#   2.  top_branch_genes_pseudotime_curves.pdf   – diagnostic plot of top genes
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(monocle3)  # ≥ 1.3.1
  library(dplyr)
  library(igraph)
  library(stringr)
  library(splines)  # for natural cubic splines in model formula
})
cds <- load_monocle_objects(directory_path=paste0(dataset,'_subset_cds.rds'))
# ==============================================================================
# 0  User‑editable parameters ----------------------------------------------------
# ==============================================================================
cds_rds_path   <- "S-MAG_cds_rooted.rds"  # path to input CDS
output_prefix  <- "MG_vs_AG"              # basename for output files

# Column names inside colData
col_gland      <- "gland"        # must contain "MG" / "AG"
col_celltype   <- "newcelltype"  # cell‑type annotation
col_stage      <- "stage"        # developmental stage (e.g. stage0 … stage4)

root_celltypes <- c("Basal", "StemCells")  # labels that define the origin
root_stage     <- "stage0"                   # earliest stage label (optional)

qvalue_cutoff  <- 0.05    # FDR threshold for significant interaction genes
n_cores        <- 1       # parallel cores for fit_models()
# ------------------------------------------------------------------------------

cat("\n===============================================\n")
cat(" Loading cell_data_set →", cds_rds_path, "\n")
cat("===============================================\n\n")
#cds <- readRDS(cds_rds_path)
cds$pseudotime <- pseudotime(cds, reduction_method = "UMAP")
# ---- sanity checks -----------------------------------------------------------
if (!"pseudotime" %in% colnames(colData(cds))) {
  stop("❌  Pseudotime not found – run order_cells() before this script.")
}
if (!all(c(col_gland, col_celltype) %in% colnames(colData(cds)))) {
  stop("❌  Required metadata columns missing in colData.")
}

# 1.2 取出 pseudotime 向量
pt <- pseudotime(cds, reduction_method = "UMAP")
# 或者直接
# pt <- colData(cds)$pseudotime

# 1.3 写成 CSV，行名是细胞条码
pt_df <- data.frame(cell = names(pt), pseudotime = pt, row.names = NULL)
write.csv(pt_df, "S_MAG_pseudotime.csv", quote = FALSE, row.names = FALSE)