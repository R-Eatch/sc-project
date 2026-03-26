#!/usr/bin/env Rscript
# ------------------------------------------------------------
# Split a Seurat object by Age_group (2 groups),
# randomly sample 10k cells from each subset (seed = 2025),
# and save two RDS files. 100% controlled via GLOBAL VARIABLES.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  suppressWarnings(library(Seurat))
})

# ===================== GLOBAL VARIABLES =====================
INPUT_RDS <- "../../data/splits_seurat/combined_seurat_full.rds"   # <<<<< set me
OUT_DIR   <- "./subset_out"               # <<<<< set me
COL_NAME  <- "Age_group"                  # metadata column with exactly 2 groups
N_TARGET  <- 20000                         # cells to sample for each group
SEED      <- 2025                          # random seed
# ============================================================

# --------------------------- helpers ------------------------
stopf    <- function(...) stop(sprintf(...), call. = FALSE)
messagef <- function(...) message(sprintf(...))

sanitize <- function(x) {
  # make safe for filenames
  x <- gsub("[\\/:*?\"<>| ]+", "_", x)
  x <- gsub("__+", "_", x)
  x <- sub("^_+", "", x)
  x <- sub("_+$", "", x)
  x
}

# --------------------------- checks -------------------------
if (!file.exists(INPUT_RDS)) stopf("Input not found: %s", INPUT_RDS)
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!is.numeric(N_TARGET) || N_TARGET <= 0) stopf("N_TARGET must be a positive integer")

set.seed(SEED)

# --------------------------- load ---------------------------
messagef("Loading Seurat object: %s", INPUT_RDS)
obj <- readRDS(INPUT_RDS)
if (!inherits(obj, "Seurat")) stopf("Input RDS is not a Seurat object (class = %s)", paste(class(obj), collapse = ","))

md <- obj@meta.data
if (!(COL_NAME %in% colnames(md))) stopf("Metadata column '%s' not found.", COL_NAME)

levels_vec <- unique(as.character(md[[COL_NAME]]))
levels_vec <- levels_vec[!is.na(levels_vec)]
if (length(levels_vec) != 2) {
  stopf("Expected exactly 2 groups in '%s', found: %s", COL_NAME, paste(levels_vec, collapse=","))
}

# --------------------------- core ---------------------------
base <- basename(INPUT_RDS)
base_noext <- sub("\\.rds$", "", base, ignore.case = TRUE)

summary_rows <- list()

for (grp in levels_vec) {
  messagef("Processing group: %s", grp)
  cells_grp <- rownames(md)[md[[COL_NAME]] == grp]
  if (length(cells_grp) == 0) {
    messagef("[WARN] No cells for group '%s' — skipping.", grp)
    next
  }

  obj_grp <- subset(obj, cells = cells_grp)
  n_avail <- ncol(obj_grp)
  n_take  <- min(N_TARGET, n_avail)
  if (n_take < N_TARGET) messagef("[INFO] Group '%s' has only %d cells; taking all.", grp, n_avail)

  sel <- sample(colnames(obj_grp), size = n_take, replace = FALSE)
  obj_grp_sub <- subset(obj_grp, cells = sel)

  grp_tag  <- sanitize(as.character(grp))
  out_file <- file.path(OUT_DIR, sprintf("%s_%s_%dk.rds", base_noext, grp_tag, as.integer(N_TARGET/1000)))
  saveRDS(obj_grp_sub, out_file)
  messagef("Saved: %s (cells=%d)", out_file, ncol(obj_grp_sub))

  summary_rows[[length(summary_rows)+1]] <- data.frame(
    group = grp,
    total_cells_in_group = n_avail,
    sampled_cells = n_take,
    output_rds = out_file,
    stringsAsFactors = FALSE
  )
}

# -------------------------- summary -------------------------
if (length(summary_rows)) {
  summary_df <- do.call(rbind, summary_rows)
  csv_path <- file.path(OUT_DIR, sprintf("%s_sampling_summary.csv", base_noext))
  utils::write.csv(summary_df, csv_path, row.names = FALSE)
  messagef("Summary written: %s", csv_path)
} else {
  messagef("No groups processed; nothing to save.")
}

messagef("Done.")

