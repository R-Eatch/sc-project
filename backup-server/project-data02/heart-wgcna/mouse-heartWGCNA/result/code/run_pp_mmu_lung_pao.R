#!/usr/bin/env Rscript
# ============================================================
# Mouse lung RNA-seq preprocessing for WGCNA (counts -> filterByExpr -> VST)
#   1) Read raw count table (featureCounts-like): Geneid + annotation cols + samples
#   2) [NEW] Keep only PAO arm samples + their Ctrl_M healthy controls
#   3) Map Ensembl gene IDs -> gene_name using a GTF
#   4) Collapse duplicate gene_name rows by sum
#   5) Low-expression filtering using edgeR::filterByExpr (counts-level)
#   6) DESeq2 VST transform (blind=TRUE) for WGCNA input
#   7) Save:
#        - <PREFIX>_counts_geneXsample_geneName_collapsed_unfiltered.csv
#        - <PREFIX>_counts_geneXsample_geneName_filtered.csv
#        - <PREFIX>_VST_geneXsample_geneName.csv
#        - <PREFIX>_datExpr_samplesXgenes_VST.csv
#        - <PREFIX>_id_mapping_report.csv
# Notes:
#   - No variance top-N filtering (per your request)
#   - Filtering is performed BEFORE VST
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# SETTINGS (EDIT ME)
# =========================
# [CHANGED] Input is CSV this time
LUNG_COUNTS <- "D:/111/mouse-heartWGCNA/data/GSE204709/GSE204709_Gene_count.csv"  # <-- your CSV/CSV.GZ counts
MOUSE_GTF   <- "D:\\111\\mouse-heartWGCNA\\data\\Mus_musculus.GRCm39.115\\Mus_musculus.GRCm39.115.gtf"

OUT_DIR <- "D:/111/mouse-heartWGCNA/data/lung_pao"
PREFIX  <- "lung_pao"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# [CHANGED] Keep only PAO arm + its healthy controls
#   Examples to KEEP:
#     Ctrl_M2_aligned.sam, Ctrl_M3_aligned.sam, ...
#     D2_PAO_M1_aligned.sam, D4_PAO_M6_aligned.sam, D14_PAO_M5_aligned.sam, ...
KEEP_SAMPLE_REGEX <- "^(Ctrl_M[0-9]+|D[0-9]+_PAO_M[0-9]+)"  # matches BEFORE suffix like _aligned.sam

# Group inference for PAO lung:
# samples whose names contain CONTROL_PATTERN are considered "control"; others are "case"
CONTROL_PATTERN <- "^Ctrl_M"   # [CHANGED] PAO healthy controls

# VST settings
VST_BLIND <- TRUE         # recommended for unsupervised use (WGCNA)

# =========================
# HELPERS
# =========================
strip_ensembl_version <- function(x) {
  ifelse(grepl("^ENS", x), sub("\\.[0-9]+$", "", x), x)
}

extract_attr <- function(attr, key) {
  attr <- as.character(attr)
  attr[is.na(attr)] <- ""
  needle <- paste0(key, ' "')
  has <- grepl(needle, attr, fixed = TRUE)
  
  out <- rep(NA_character_, length(attr))
  if (!any(has)) return(out)
  
  pat_full <- paste0('.*', key, ' "([^"]+)".*')
  out[has] <- sub(pat_full, "\\1", attr[has], perl = TRUE)
  
  bad <- has & (out[has] == attr[has])
  out[bad] <- NA_character_
  out
}

infer_group_from_colnames <- function(sample_names, control_pattern = "TC") {
  is_control <- grepl(control_pattern, sample_names, ignore.case = TRUE)
  factor(ifelse(is_control, "control", "case"), levels = c("control", "case"))
}

read_gtf_gene_map <- function(gtf_path) {
  message("[GTF] Reading gene map from: ", gtf_path)
  if (!file.exists(gtf_path)) stop("GTF not found: ", gtf_path)
  
  find_gtf_skip <- function(path, chunk = 20000L, max_lines = 500000L) {
    con <- file(path, open = "r")
    on.exit(close(con), add = TRUE)
    skipped <- 0L
    repeat {
      lines <- readLines(con, n = chunk, warn = FALSE)
      if (length(lines) == 0L) break
      for (ln in lines) {
        if (!nzchar(ln)) { skipped <- skipped + 1L; next }
        if (startsWith(ln, "#")) {
          skipped <- skipped + 1L
        } else {
          return(skipped)
        }
        if (skipped >= max_lines) break
      }
      if (skipped >= max_lines) break
    }
    stop("Could not find first non-comment line in GTF within ", max_lines, " lines: ", path)
  }
  
  skip_n <- find_gtf_skip(gtf_path)
  message("[GTF] Skipping ", skip_n, " comment/header lines starting with '#'")
  
  gtf <- fread(
    gtf_path,
    sep = "\t",
    header = FALSE,
    quote = "",
    data.table = TRUE,
    fill = TRUE,
    skip = skip_n
  )
  if (ncol(gtf) < 9) stop("GTF appears malformed (<9 columns): ", gtf_path)
  setnames(gtf, paste0("V", seq_len(ncol(gtf))))
  
  gene_rows <- gtf[V3 == "gene"]
  if (nrow(gene_rows) == 0) {
    message("[GTF] No 'gene' rows found; falling back to all rows.")
    gene_rows <- gtf
  }
  
  attr <- gene_rows$V9
  gene_id   <- strip_ensembl_version(extract_attr(attr, "gene_id"))
  gene_name <- extract_attr(attr, "gene_name")
  
  map <- data.table(gene_id = gene_id, gene_name = gene_name)
  map <- map[!is.na(gene_id) & gene_id != ""]
  setorder(map, gene_id)
  
  # Deduplicate gene_id -> first non-empty gene_name
  map <- map[, .(gene_name = {
    x <- gene_name[!is.na(gene_name) & gene_name != ""]
    if (length(x) == 0) NA_character_ else x[1]
  }), by = gene_id]
  
  message("[GTF] Mapped gene_id entries: ", nrow(map))
  map
}

# =========================
# RUN
# =========================
message("========== [1] Read raw counts ==========")
# fread auto-detects CSV/TSV; works for .csv/.csv.gz
dt <- fread(LUNG_COUNTS)
if (nrow(dt) == 0) stop("Empty input: ", LUNG_COUNTS)

# Detect gene id column
if ("Geneid" %in% names(dt)) {
  gene_id_raw <- dt[["Geneid"]]
  annot_cols <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
} else if ("gene_id" %in% names(dt)) {
  gene_id_raw <- dt[["gene_id"]]
  annot_cols <- "gene_id"
} else {
  gene_id_raw <- dt[[1]]
  annot_cols <- names(dt)[1]
  message("[WARN] Could not find 'Geneid' or 'gene_id'; using first column as gene id: ", annot_cols)
}

gene_id_raw <- as.character(gene_id_raw)
gene_id_nv  <- strip_ensembl_version(gene_id_raw)

count_cols <- setdiff(names(dt), annot_cols)
if (length(count_cols) == 0) stop("No sample count columns detected. Columns: ", paste(names(dt), collapse = ", "))

# ========== [NEW] Keep only PAO + Ctrl_M samples ==========
keep_cols <- count_cols[grepl(KEEP_SAMPLE_REGEX, count_cols, perl = TRUE)]
if (length(keep_cols) == 0) {
  stop("No columns matched KEEP_SAMPLE_REGEX=", KEEP_SAMPLE_REGEX,
       "\nAvailable sample columns preview: ", paste(head(count_cols, 30), collapse = ", "))
}

# Optional: keep deterministic order
keep_cols <- keep_cols[order(keep_cols)]

message("========== [1b] Subset samples (PAO + Ctrl_M) ==========")
message("Keeping samples: ", length(keep_cols), " / ", length(count_cols))
message("Kept sample columns preview:")
print(head(keep_cols, 30))

count_mat <- as.matrix(dt[, ..keep_cols])
storage.mode(count_mat) <- "double"
rownames(count_mat) <- gene_id_nv

message("Counts dim (genes x samples): ", paste(dim(count_mat), collapse = " x "))

message("========== [2] Map Ensembl gene_id -> gene_name (GTF) ==========")
map <- read_gtf_gene_map(MOUSE_GTF)
m <- map$gene_name; names(m) <- map$gene_id

gene_name <- unname(m[gene_id_nv])
n_total <- length(gene_id_nv)
n_mapped <- sum(!is.na(gene_name) & gene_name != "")
message("Mapped gene_name: ", n_mapped, " / ", n_total)

report <- data.table(
  gene_id_input = gene_id_raw,
  gene_id_noversion = gene_id_nv,
  gene_name = gene_name
)

use_name <- gene_name
use_name[is.na(use_name) | use_name == ""] <- gene_id_nv[is.na(use_name) | use_name == ""]

message("========== [3] Collapse duplicated gene names by sum ==========")
DT <- as.data.table(count_mat)
DT[, gene := use_name]
DT2 <- DT[, lapply(.SD, sum), by = gene, .SDcols = colnames(count_mat)]

count_mat2 <- as.matrix(DT2[, -1, with = FALSE])
rownames(count_mat2) <- DT2$gene
storage.mode(count_mat2) <- "double"
message("After collapse dim (genes x samples): ", paste(dim(count_mat2), collapse = " x "))

# Save collapsed-unfiltered counts for traceability
out_counts_unf <- file.path(OUT_DIR, paste0(PREFIX, "_counts_geneXsample_geneName_collapsed_unfiltered.csv"))
fwrite(as.data.table(count_mat2, keep.rownames = "gene"), out_counts_unf)

message("========== [4] Low-expression filtering: edgeR::filterByExpr ==========")
if (!requireNamespace("edgeR", quietly = TRUE)) {
  stop("edgeR is required for filterByExpr but is not installed. Please install edgeR.")
}

group <- infer_group_from_colnames(colnames(count_mat2), CONTROL_PATTERN)
message("Group counts: ", paste(names(table(group)), as.integer(table(group)), sep = "=", collapse = ", "))

# Drop all-zero genes first (always safe)
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

out_counts_filt <- file.path(OUT_DIR, paste0(PREFIX, "_counts_geneXsample_geneName_filtered.csv"))
fwrite(as.data.table(count_filt, keep.rownames = "gene"), out_counts_filt)

message("========== [5] DESeq2 VST ==========")
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  stop("DESeq2 is required for VST but is not installed. Please install DESeq2.")
}

meta <- data.frame(row.names = colnames(count_filt), group = group)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_filt),
  colData   = meta,
  design    = ~ group
)

dds <- DESeq2::estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind = VST_BLIND)

if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  stop("SummarizedExperiment is required (it is a DESeq2 dependency) but not installed.")
}
vst_mat <- SummarizedExperiment::assay(vsd)
message("VST dim (genes x samples): ", paste(dim(vst_mat), collapse = " x "))
message("Any NA in VST: ", anyNA(vst_mat))
message("Any non-finite in VST: ", any(!is.finite(vst_mat)))
message("Range of VST: ", paste(range(vst_mat, finite = TRUE), collapse = " ~ "))

message("========== [6] Build WGCNA datExpr (samples x genes) ==========")
datExpr <- t(vst_mat)
message("datExpr dim (samples x genes): ", paste(dim(datExpr), collapse = " x "))

message("========== [7] Save outputs ==========")
out_vst     <- file.path(OUT_DIR, paste0(PREFIX, "_VST_geneXsample_geneName.csv"))
out_datExpr <- file.path(OUT_DIR, paste0(PREFIX, "_datExpr_samplesXgenes_VST.csv"))
out_report  <- file.path(OUT_DIR, paste0(PREFIX, "_id_mapping_report.csv"))

fwrite(as.data.table(vst_mat, keep.rownames = "gene"), out_vst)
fwrite(as.data.table(datExpr, keep.rownames = "sample"), out_datExpr)
fwrite(report, out_report)

message("Done.")
message("  Collapsed counts (unfiltered): ", out_counts_unf)
message("  Filtered counts (filterByExpr): ", out_counts_filt)
message("  VST matrix: ", out_vst)
message("  WGCNA datExpr (VST): ", out_datExpr)
message("  Mapping report: ", out_report)
