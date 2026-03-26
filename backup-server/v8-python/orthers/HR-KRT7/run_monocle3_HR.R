library(Seurat)
library(monocle3)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)

#### 全局变量 ####
DATASETS <- c("M-MG","S-MG")
cores <- 1
DRAW_PLOTS <- TRUE
do_preprocess <-FALSE
# ---- Paths aligned with previous step ----
in_dir  <- file.path("data")
out_dir <- file.path("result", "monocle3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

process_one_dataset <- function(dataset, cores = 1) {
  in_rds <- file.path(in_dir, paste0(dataset, "_filtered_reumap.rds"))
  out_cds_rds <- file.path(out_dir, paste0(dataset, "_cds_for_subset_HR.rds"))

  if (!file.exists(in_rds)) {
    message("[SKIP] missing input RDS: ", in_rds)
    return(invisible(NULL))
  }
  if (file.exists(out_cds_rds)) {
    message("[SKIP] exists: ", out_cds_rds)
    return(invisible(NULL))
  }

  message("\n", paste(rep("=", 80), collapse=""))
  message("[RUN] ", dataset)

  ea <- readRDS(in_rds)

  data <- GetAssayData(ea, assay = "RNA", slot = "counts")
  cell_metadata <- ea@meta.data
  gene_annotation <- data.frame(gene_short_name = rownames(data))
  rownames(gene_annotation) <- rownames(data)

  cds <- new_cell_data_set(
    data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )

  rm(data, cell_metadata, gene_annotation)

  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, cores = cores)

  # 覆盖 monocle 的 UMAP，用 Seurat 的 UMAP 坐标（保持和上一步一致的嵌入）
  int.embed <- Embeddings(ea, reduction = "umap")
  cds.embed <- cds@int_colData$reducedDims$UMAP
  int.embed <- int.embed[rownames(cds.embed), ]
  cds@int_colData$reducedDims$UMAP <- int.embed

  # 可选：画图检查（不保存）
  plot_cells(cds, reduction_method = "UMAP", color_cells_by = "newcelltype")

  saveRDS(cds, file = out_cds_rds)
  message("[OK] saved: ", out_cds_rds)

  invisible(TRUE)
}

# ---- Batch run ----
if(do_preprocess){
  for (ds in DATASETS) {
  process_one_dataset(ds, cores = cores)
  }
}

message("\nDone.")

# aligned paths
in_dir  <- file.path("result", "monocle3")
out_dir <- file.path("result", "monocle3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# =========================
# Functions
# =========================
#dataset <- 'M-MG'
if(FALSE){
  cds_subset <- load_monocle_objects(directory_path=paste0('D:/111/',dataset,'_cds_for_subset_HR.rds'))
  cds_subset <- cluster_cells(cds_subset,random_seed = 2024)
  cds_subset <- learn_graph(cds_subset,use_partition = FALSE)
  cds_subset <- order_cells(cds_subset)   ####手动操作####
  saveRDS(cds_subset, file=paste0('D:/111/',dataset,'_subset_cds_HR.rds'))
}

process_one_dataset <- function(dataset) {

  in_cds_rds  <- file.path(in_dir, paste0(dataset, "_subset_cds_HR.rds"))


  out_pt_pdf  <- file.path(out_dir, paste0(dataset, "_m3_pseudotime.pdf"))
  out_ct_pdf  <- file.path(out_dir, paste0(dataset, "_m3_newcelltype.pdf"))



  cds <- readRDS(in_cds_rds)


  # 4) plots
  p_ct <- plot_cells(
    cds,
    color_cells_by = "newcelltype"
                      #"stage"
,
    label_groups_by_cluster = TRUE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 1
  )
  ggsave(out_ct_pdf, p_ct, width = 3, height = 3)
  message("[SAVE] ", out_ct_pdf)

    p_pt <- plot_cells(
      cds,
      color_cells_by = "pseudotime",
      label_cell_groups = TRUE,
      label_leaves = TRUE,
      label_branch_points = TRUE,
      graph_label_size = 1
    )
    ggsave(out_pt_pdf, p_pt, width = 3, height = 3)
    message("[SAVE] ", out_pt_pdf)

  invisible(TRUE)
}

# =========================
# Batch run
# =========================
if(DRAW_PLOTS){
  for (ds in DATASETS) {
  process_one_dataset(ds)
  }
    message("\nDone.")
}else{
  message("\nDo not draw plot.")
}

