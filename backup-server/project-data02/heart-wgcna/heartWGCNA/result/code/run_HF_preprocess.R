#!/usr/bin/env Rscript
# ============================================================
# Preprocess multi-chamber HF expression matrix -> WGCNA input
# Output as CSV (NOT RDS), per your request.
#
# Steps:
#   1) Keep ALL samples (no chamber filtering)
#   2) Remove genes that are all-NA across samples
#   3) Fill remaining NA with 0
#   4) Output WGCNA input matrix: samples x genes (CSV)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# User settings
# -----------------------------
infile  <- "D:\\111\\heartWGCNA\\data\\GSE161472_FPKM\\GSE161472_FPKM.txt"   # <-- change if needed
outdir  <- "D:\\111\\heartWGCNA\\data\\GSE161472_FPKM\\preprocessed"          # <-- change if needed
prefix  <- "GSE161472_ALL_preprocessed"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) Read matrix (quoted, NA present)
# -----------------------------
message("[1/5] Reading input matrix...")
dt <- fread(
  infile,
  sep = "\t",
  header = TRUE,
  quote = "\"",
  na.strings = c("NA", ""),
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(dt) < 2) stop("Input file seems to have <2 columns. Check delimiter/encoding.")

gene_id <- dt[[1]]
expr0   <- as.matrix(dt[, -1, drop = FALSE])
rownames(expr0) <- gene_id

message("  Raw dim (genes x samples): ", paste(dim(expr0), collapse = " x "))

# -----------------------------
# 2) Keep ALL columns (no filtering)
# -----------------------------
message("[2/5] Keeping ALL samples (no chamber filtering)...")
expr_all <- expr0

message("  ALL dim (genes x samples): ", paste(dim(expr_all), collapse = " x "))

# -----------------------------
# 3) Remove genes that are all NA across samples
# -----------------------------
message("[3/5] Removing all-NA genes...")
all_na <- rowSums(is.na(expr_all)) == ncol(expr_all)
message("  All-NA genes to drop: ", sum(all_na), " / ", nrow(expr_all))
expr_all <- expr_all[!all_na, , drop = FALSE]

# -----------------------------
# 4) Fill remaining NA with 0
# -----------------------------
message("[4/5] Filling remaining NA with 0...")
na_count <- sum(is.na(expr_all))
message("  NA to fill: ", na_count)
if (na_count > 0) expr_all[is.na(expr_all)] <- 0

storage.mode(expr_all) <- "double"

# -----------------------------
# 5) Build WGCNA input matrix (samples x genes) and save as CSV
# -----------------------------
message("[5/5] Building WGCNA input (samples x genes) and saving CSV...")
datExpr <- t(expr_all)  # samples x genes

message("  datExpr dim (samples x genes): ", paste(dim(datExpr), collapse = " x "))
message("  Any NA left? ", any(is.na(datExpr)))
message("  Any non-finite? ", any(!is.finite(datExpr)))

# Write CSV (rowname as first column)
out_csv <- file.path(outdir, paste0(prefix, "_datExpr_samplesXgenes.csv"))
dat_out <- data.table::as.data.table(datExpr, keep.rownames = "sample")
fwrite(dat_out, out_csv, sep = ",", quote = TRUE)

# (Optional) also save genes x samples for debugging (CSV)
out_csv2 <- file.path(outdir, paste0(prefix, "_expr_genesXsamples.csv"))
expr_out <- data.table::as.data.table(expr_all, keep.rownames = "gene")
fwrite(expr_out, out_csv2, sep = ",", quote = TRUE)

message("Done.")
message("WGCNA input CSV: ", out_csv)
message("Reference CSV (genes x samples): ", out_csv2)
