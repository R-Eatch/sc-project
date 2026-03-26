#!/usr/bin/env Rscript
# ============================================================
# Mouse HEART RNA-seq preprocessing for WGCNA (counts -> filterByExpr -> VST)
# Input format (already gene symbols):
#   Genesymbol  sample1 sample2 ...
#   Xkr4        1       2       ...
# Output:
#   1) <PREFIX>_counts_geneXsample_geneSymbol_collapsed_unfiltered.csv
#   2) <PREFIX>_counts_geneXsample_geneSymbol_filtered.csv
#   3) <PREFIX>_VST_geneXsample_geneSymbol.csv
#   4) <PREFIX>_datExpr_samplesXgenes_VST.csv   (WGCNA input)
# Notes:
#   - No gene_id->gene_name mapping needed (already gene symbols)
#   - Collapses duplicated gene symbols by sum (recommended)
#   - Uses lung-style preprocessing: filterByExpr (edgeR) then VST (DESeq2)
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# SETTINGS (EDIT ME)
# =========================
HEART_COUNTS <- "D:/111/mouse-heartWGCNA/data/GSE236374_MI_FeatureCounts_matrix/GSE236374_MI_FeatureCounts_matrix.txt"  # <-- change
OUT_DIR      <- "D:/111/mouse-heartWGCNA/data/heart"                               # <-- change
PREFIX       <- "heart"                                                            # <-- change

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Group inference for heart (per your earlier convention):
# samples whose names contain CONTROL_PATTERN are considered "control"; others are "case"
CONTROL_PATTERN <- "sham"  # change if needed

# VST settings
VST_BLIND <- TRUE          # recommended for unsupervised use (WGCNA)

# =========================
# HELPERS
# =========================
clean_gene_symbol <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  x
}

infer_group_from_colnames <- function(sample_names, control_pattern = "sham") {
  is_control <- grepl(control_pattern, sample_names, ignore.case = TRUE)
  factor(ifelse(is_control, "control", "case"), levels = c("control", "case"))
}

# =========================
# RUN
# =========================
message("========== [1] Read HEART counts (gene symbols) ==========")
dt <- fread(HEART_COUNTS)
if (nrow(dt) == 0) stop("Empty input: ", HEART_COUNTS)

# Detect gene symbol column
sym_col <- NULL
if ("Genesymbol" %in% names(dt)) sym_col <- "Genesymbol"
if (is.null(sym_col) && "gene" %in% names(dt)) sym_col <- "gene"
if (is.null(sym_col) && "Gene" %in% names(dt)) sym_col <- "Gene"
if (is.null(sym_col) && "symbol" %in% names(dt)) sym_col <- "symbol"
if (is.null(sym_col)) sym_col <- names(dt)[1]
message("Gene symbol column: ", sym_col)

symbols <- clean_gene_symbol(dt[[sym_col]])
if (all(is.na(symbols))) stop("All gene symbols are NA/blank after cleaning.")

# Sample columns = all except sym_col
sample_cols <- setdiff(names(dt), sym_col)
if (length(sample_cols) == 0) stop("No sample columns found. Columns: ", paste(names(dt), collapse = ", "))

# Coerce counts to numeric safely
count_mat <- as.matrix(dt[, ..sample_cols])
suppressWarnings(storage.mode(count_mat) <- "double")

nonfinite_n <- sum(!is.finite(count_mat))
if (nonfinite_n > 0) {
  message("[WARN] Non-finite values detected in count matrix: ", nonfinite_n, ". Converting them to NA then 0.")
  count_mat[!is.finite(count_mat)] <- NA_real_
}
count_mat[is.na(count_mat)] <- 0

rownames(count_mat) <- symbols

message("Counts dim (genes x samples): ", paste(dim(count_mat), collapse = " x "))
message("Sample preview: ")
print(head(sample_cols, 30))

message("========== [2] Collapse duplicated gene symbols by sum ==========")
DT <- as.data.table(count_mat)
DT[, gene := rownames(count_mat)]
DT2 <- DT[, lapply(.SD, sum), by = gene, .SDcols = colnames(count_mat)]

count_mat2 <- as.matrix(DT2[, -1, with = FALSE])
rownames(count_mat2) <- DT2$gene
storage.mode(count_mat2) <- "double"

message("After collapse dim (genes x samples): ", paste(dim(count_mat2), collapse = " x "))

# Quick QC: library sizes
libsize <- colSums(count_mat2)
message("Library size summary (raw counts):")
print(summary(libsize))

# Save collapsed-unfiltered counts for traceability
out_counts_unf <- file.path(OUT_DIR, paste0(PREFIX, "_counts_geneXsample_geneSymbol_collapsed_unfiltered.csv"))
fwrite(as.data.table(count_mat2, keep.rownames = "gene"), out_counts_unf)

message("========== [3] Low-expression filtering: edgeR::filterByExpr ==========")
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is required for filterByExpr but is not installed. Please install edgeR.")
}

group <- infer_group_from_colnames(colnames(count_mat2), CONTROL_PATTERN)
message("Group counts: ", paste(names(table(group)), as.integer(table(group)), sep = "=", collapse = ", "))

# Drop all-zero genes first
nonzero <- rowSums(count_mat2) > 0
if (any(!nonzero)) {
  message("Dropping all-zero genes: ", sum(!nonzero))
  count_mat2 <- count_mat2[nonzero, , drop = FALSE]
}

y0 <- edgeR::DGEList(counts = round(count_mat2), group = group)
keep <- edgeR::filterByExpr(y0, group = group)
message("filterByExpr kept genes: ", sum(keep), " / ", length(keep))

count_filt <- y0$counts[keep, , drop = FALSE]
message("After filterByExpr dim (genes x samples): ", paste(dim(count_filt), collapse = " x "))

out_counts_filt <- file.path(OUT_DIR, paste0(PREFIX, "_counts_geneXsample_geneSymbol_filtered.csv"))
fwrite(as.data.table(count_filt, keep.rownames = "gene"), out_counts_filt)

message("========== [4] DESeq2 VST ==========")
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required for VST but is not installed. Please install DESeq2.")
}
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("SummarizedExperiment is required (DESeq2 dependency) but is not installed.")
}

meta <- data.frame(row.names = colnames(count_filt), group = group)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_filt),
  colData   = meta,
  design    = ~ group
)

dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = VST_BLIND)

vst_mat <- SummarizedExperiment::assay(vsd)  # genes x samples
message("VST dim (genes x samples): ", paste(dim(vst_mat), collapse = " x "))
message("Any NA in VST: ", anyNA(vst_mat))
message("Any non-finite in VST: ", any(!is.finite(vst_mat)))
message("Range of VST: ", paste(range(vst_mat, finite = TRUE), collapse = " ~ "))

message("========== [5] Build WGCNA datExpr (samples x genes) ==========")
datExpr <- t(vst_mat)
message("datExpr dim (samples x genes): ", paste(dim(datExpr), collapse = " x "))

message("========== [6] Save outputs ==========")
out_vst     <- file.path(OUT_DIR, paste0(PREFIX, "_VST_geneXsample_geneSymbol.csv"))
out_datExpr <- file.path(OUT_DIR, paste0(PREFIX, "_datExpr_samplesXgenes_VST.csv"))

fwrite(as.data.table(vst_mat, keep.rownames = "gene"), out_vst)
fwrite(as.data.table(datExpr, keep.rownames = "sample"), out_datExpr)

message("Done.")
message("  Collapsed counts (unfiltered): ", out_counts_unf)
message("  Filtered counts (filterByExpr): ", out_counts_filt)
message("  VST matrix: ", out_vst)
message("  WGCNA datExpr: ", out_datExpr)
