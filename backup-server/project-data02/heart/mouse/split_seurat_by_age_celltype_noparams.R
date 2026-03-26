#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# split_seurat_by_age_celltype_noparams.R
# Workflow:
# 1) Convert input .h5ad (subset) -> Seurat (via sceasy)
# 2) If assay@data appears empty (i.e., only raw counts), run NormalizeData
# 3) HVGs -> ScaleData -> PCA -> FindNeighbors -> FindClusters -> UMAP
# 4) Save two UMAP figures: colored by Main_cell_type and by Age_group
# 5) Split into 6 RDS by Age_group × Main_cell_type
#
# All paths/columns are hard-coded below.

suppressPackageStartupMessages({
  library(reticulate)
  library(sceasy)
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

# ---------------- Global config ----------------
#INPUT_H5AD <- "../data/GSE247719_PanSci_03_Heart_adata.sub_Age03_23_Fibro_VentCM_EndoEndo.h5ad"
INPUT_H5AD <-"../data/GSE247719_PanSci_03_Heart_adata.sub_Age03_23_Main_cells.h5ad"

OUTDIR     <- "../data/splits_seurat"
USE_RAW    <- FALSE     # If TRUE, copy .raw to X before converting (via scanpy/anndata)
PYTHON_BIN <- NULL      # e.g., "/path/to/conda/envs/yourenv/bin/python" or NULL to auto
AGE_COL    <- "Age_group"
TYPE_COL   <- "Main_cell_type"
AGE_LEVELS <- c("03_months", "23_months")
#TYPE_LEVELS<- c("Fibroblasts", "Ventricular cardiomyocytes", "Endocardial endothelial cells")
TYPE_LEVELS<- c(
    "Vascular endothelial cells",
    "Fibroblasts",
    "Ventricular cardiomyocytes",
    "Myeloid cells",
    "Mural cells",
    "Endocardial endothelial cells",
    "Lymphoid cells_T cells",
    "Lymphatic endothelial cells",
    "Atrial cardiomyocytes",
    "Lymphoid cells_B cells",
    "Mesothelial cells",
    "Neural cells",
    "Adipocytes"
)
# Dimensionality & graph settings
N_PCS      <- 20
K_PARAM    <- 30
RESOLUTION <- 0.4
UMAP_MIN_DIST <- 0.3
UMAP_N_NEIGHBORS <- 30
FIG_WIDTH  <- 8
FIG_HEIGHT <- 6
FIG_DPI    <- 300
SEED       <- 12345
# -----------------------------------------------

set.seed(SEED)

if (!is.null(PYTHON_BIN)) {
  reticulate::use_python(PYTHON_BIN, required = TRUE)
  message(sprintf("Using Python: %s", PYTHON_BIN))
}

if (!file.exists(INPUT_H5AD)) {
  stop(sprintf("Input file not found: %s", INPUT_H5AD))
}
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Optionally generate a temp h5ad with raw->X (if USE_RAW)
tmpfile <- INPUT_H5AD
if (USE_RAW) {
  message("USE_RAW=TRUE: reading via scanpy and copying .raw to X ...")
  scanpy <- reticulate::import("scanpy", convert = FALSE)
  ad <- scanpy$read_h5ad(INPUT_H5AD)
  if (is.null(reticulate::py_get_attr(ad, "raw"))) stop("Requested USE_RAW but adata.raw is NULL / missing.")
  ad$X <- ad$raw$X$copy()
  ad$var <- ad$raw$var$copy()
  tmpfile <- tempfile(pattern="h5ad_raw_", fileext = ".h5ad")
  ad$write_h5ad(tmpfile, compression = "gzip")
  message(sprintf("Temporary raw->X h5ad: %s", tmpfile))
}

message(sprintf("Converting to Seurat: %s", INPUT_H5AD))
obj <- sceasy::convertFormat(tmpfile, from = "anndata", to = "seurat")

# Ensure metadata presence
meta <- obj@meta.data
need_cols <- c(AGE_COL, TYPE_COL)
if (!all(need_cols %in% colnames(meta))) {
  stop(sprintf("Missing required metadata columns: %s ; have: %s",
               paste(setdiff(need_cols, colnames(meta)), collapse=", "),
               paste(colnames(meta), collapse=", ")))
}
# Coerce to factor for nice plotting
for (cn in need_cols) {
  if (!is.factor(meta[[cn]])) meta[[cn]] <- factor(meta[[cn]])
}
obj@meta.data <- meta

# Determine if normalization needed
assay <- DefaultAssay(obj)
data_mat   <- GetAssayData(obj, assay = assay, slot = "data")
nnz_data <- if (is(data_mat, "dgCMatrix")) Matrix::nnzero(data_mat) else sum(abs(as.matrix(data_mat)) != 0)
need_norm <- (ncol(data_mat) == 0) || (nnz_data == 0)

message(sprintf("Assay: %s | data nnz = %s | need_norm = %s", assay, format(nnz_data, big.mark=","), need_norm))
need_norm=TRUE
if (need_norm) {
  message("Running NormalizeData (LogNormalize, scale.factor=1e4) ...")
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
} else {
  message("Data appears normalized; skipping NormalizeData.")
}

# Feature selection, scaling, PCA
message("Finding variable features (vst, 3000) ...")
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)

message("Scaling data (variable features) ...")
obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)

message(sprintf("Running PCA (npcs=%d) ...", N_PCS))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = N_PCS, verbose = FALSE)

# Neighbors & clustering
message(sprintf("Finding neighbors (dims 1:%d, k.param=%d) ...", N_PCS, K_PARAM))
obj <- FindNeighbors(obj, dims = 1:N_PCS, k.param = K_PARAM, verbose = FALSE)

message(sprintf("Finding clusters (resolution=%.2f) ...", RESOLUTION))
obj <- FindClusters(obj, resolution = RESOLUTION, algorithm = 1, verbose = FALSE)  # Louvain

# UMAP
message(sprintf("Running UMAP (dims=1:%d, n.neighbors=%d, min.dist=%.2f) ...",
                N_PCS, UMAP_N_NEIGHBORS, UMAP_MIN_DIST))
obj <- RunUMAP(obj, dims = 1:N_PCS, n.neighbors = UMAP_N_NEIGHBORS,
               min.dist = UMAP_MIN_DIST, verbose = FALSE, seed.use = SEED)

# Save UMAP figures
p1 <- DimPlot(obj, reduction = "umap", group.by = TYPE_COL, label = TRUE, repel = TRUE) +
  ggtitle(sprintf("UMAP by %s", TYPE_COL))
p2 <- DimPlot(obj, reduction = "umap", group.by = AGE_COL,  label = TRUE, repel = TRUE) +
  ggtitle(sprintf("UMAP by %s", AGE_COL))

fn_by <- function(x) gsub("\\s+", "_", tolower(x))
fig1 <- file.path(OUTDIR, sprintf("umap_by_%s.png", fn_by(TYPE_COL)))
fig2 <- file.path(OUTDIR, sprintf("umap_by_%s.png", fn_by(AGE_COL)))
ggsave(fig1, plot = p1, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
ggsave(fig2, plot = p2, width = FIG_WIDTH, height = FIG_HEIGHT, dpi = FIG_DPI)
message(sprintf("[OK] Saved figures:\n- %s\n- %s", fig1, fig2))

# Helper to make clean file names
clean <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

# Split into 6 subsets and write RDS
# dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
# n_written <- 0L
# for (a in AGE_LEVELS) {
#   for (t in TYPE_LEVELS) {
#     idx <- which(obj@meta.data[[AGE_COL]] == a & obj@meta.data[[TYPE_COL]] == t)
#     if (length(idx) == 0) {
#       message(sprintf("[SKIP] No cells for %s + %s", a, t))
#       next
#     }
#     sub_obj <- obj[, idx, drop = FALSE]
#     sub_obj@meta.data[[AGE_COL]]  <- droplevels(sub_obj@meta.data[[AGE_COL]])
#     sub_obj@meta.data[[TYPE_COL]] <- droplevels(sub_obj@meta.data[[TYPE_COL]])
#     fname <- sprintf("seurat_%s_%s.rds", clean(a), clean(t))
#     fpath <- file.path(OUTDIR, fname)
#     saveRDS(sub_obj, file = fpath)
#     message(sprintf("[OK] Wrote: %s (%d cells)", fpath, ncol(sub_obj)))
#     n_written <- n_written + 1L
#   }
# }
# message(sprintf("Finished. Total subsets written: %d", n_written))

# Save combined Seurat object as well
saveRDS(obj, file = file.path(OUTDIR, "combined_seurat_full.rds"))
message("[OK] Saved combined object.")
