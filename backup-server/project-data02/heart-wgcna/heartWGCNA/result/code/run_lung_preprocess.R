suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# SETTINGS
# =========================
IN_DIR  <- "D:/111/heartWGCNA/data/GSE262433_RAW"
OUT_DIR <- file.path(IN_DIR, "step1_merge_counts")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# STAR ReadsPerGene columns:
# 1=gene_id, 2=unstranded, 3=stranded-forward, 4=stranded-reverse
COUNT_COL <- 2

# if you want to temporarily exclude some samples by filename patterns
EXCLUDE_PATTERNS <- character(0)   # e.g. c("PL-DMEM")

# =========================
# HELPERS
# =========================
get_sample_name <- function(fp) {
  x <- basename(fp)
  x <- sub("\\.gz$", "", x)
  x <- sub("ReadsPerGene\\.out\\.tab$", "", x)
  x <- sub("_+$", "", x)
  x
}

strip_ensg_version <- function(x) {
  # ENSG00000123456.7 -> ENSG00000123456
  sub("^(ENSG[0-9]+)\\.[0-9]+$", "\\1", x)
}

read_star_counts <- function(fp, count_col = 2) {
  dt <- fread(fp, header = FALSE, sep = "\t", data.table = FALSE)
  if (ncol(dt) < 4) stop("Unexpected STAR counts format (<4 cols): ", fp)
  
  # drop STAR summary rows
  dt <- dt[!grepl("^N_", dt[[1]]), , drop = FALSE]
  
  out <- dt[, c(1, count_col), drop = FALSE]
  colnames(out) <- c("gene_id", "count")
  out
}

# =========================
# RUN
# =========================
message("========== [1] Find files ==========")
files <- list.files(IN_DIR, pattern = "ReadsPerGene\\.out\\.tab(\\.gz)?$", full.names = TRUE)
message("Found: ", length(files))
if (length(files) == 0) stop("No ReadsPerGene.out.tab(.gz) found in: ", IN_DIR)

if (length(EXCLUDE_PATTERNS) > 0) {
  drop <- Reduce(`|`, lapply(EXCLUDE_PATTERNS, function(p) grepl(p, basename(files))))
  message("Drop by pattern: ", sum(drop))
  files <- files[!drop]
}
message("Kept files: ", length(files))

sample_names <- vapply(files, get_sample_name, character(1))
if (any(duplicated(sample_names))) {
  stop("Duplicated sample names parsed: ", paste(unique(sample_names[duplicated(sample_names)]), collapse = ", "))
}
message("Sample name preview:")
print(head(sample_names, 12))

message("========== [2] Read + merge ==========")
lst <- vector("list", length(files))
names(lst) <- sample_names

for (i in seq_along(files)) {
  fp <- files[i]
  s  <- sample_names[i]
  message("Reading: ", s)
  one <- read_star_counts(fp, COUNT_COL)
  one[[s]] <- one$count
  one$count <- NULL
  lst[[s]] <- one
}

merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), lst)
mat <- as.matrix(merged[, -1, drop = FALSE])
rownames(mat) <- merged$gene_id
mat[is.na(mat)] <- 0
storage.mode(mat) <- "double"

message("Merged dim (rows x cols) = genes x samples: ", paste(dim(mat), collapse = " x "))

# Optional: keep host ENSG only
keep_ensg <- grepl("^ENSG", rownames(mat))
message("Rows starting with ENSG: ", sum(keep_ensg), " / ", nrow(mat))

# Strip ENSG version and collapse duplicates (recommended even at step1)
ensg_nv <- strip_ensg_version(rownames(mat))
if (any(duplicated(ensg_nv))) {
  message("Duplicated ENSG after stripping version: collapsing by sum.")
  dt <- as.data.table(mat)
  dt[, gene_id := ensg_nv]
  dt2 <- dt[, lapply(.SD, sum), by = gene_id, .SDcols = colnames(mat)]
  mat2 <- as.matrix(dt2[, -1, drop = FALSE])
  rownames(mat2) <- dt2$gene_id
  storage.mode(mat2) <- "double"
} else {
  mat2 <- mat
  rownames(mat2) <- ensg_nv
}

message("After strip+collapse dim: ", paste(dim(mat2), collapse = " x "))
message("Preview (first 5 genes x first 6 samples):")
print(mat2[1:min(5,nrow(mat2)), 1:min(6,ncol(mat2)), drop = FALSE])

# quick QC
libsize <- colSums(mat2)
message("Library size summary:")
print(summary(libsize))

# Save
out_counts <- file.path(OUT_DIR, "GSE262433_counts_geneXsample_ENSG_noversion.csv")
fwrite(as.data.table(mat2, keep.rownames = "gene_id"), out_counts, sep = ",", quote = TRUE)
message("Saved merged counts: ", out_counts)


############################################################################################################


suppressPackageStartupMessages({
  library(data.table)
})

# =========================
# SETTINGS
# =========================
STEP1_COUNTS <- "D:/111/heartWGCNA/data/GSE262433_RAW/step1_merge_counts/GSE262433_counts_geneXsample_ENSG_noversion.csv"
OUT_DIR      <- dirname(STEP1_COUNTS)
PREFIX       <- "GSE262433"

# Recommended for WGCNA: drop very lowly expressed genes (optional but strongly suggested)
MIN_CPM     <- 1
MIN_SAMPLES <- 3

# Optional: keep top variable genes to stabilize WGCNA (esp. small sample size)
TOP_VAR_GENES <- NA  # set NA to disable

# Prefer edgeR for TMM+logCPM; fallback to simple log2(CPM+1)
USE_EDGE_R <- TRUE

# =========================
# RUN
# =========================
message("========== [1] Read merged counts ==========")
dt <- fread(STEP1_COUNTS)
stopifnot("gene_id" %in% colnames(dt))

gene_id <- dt$gene_id
count_mat <- as.matrix(dt[, -1, with = FALSE])
rownames(count_mat) <- gene_id
storage.mode(count_mat) <- "double"

message("Counts dim (genes x samples): ", paste(dim(count_mat), collapse = " x "))
message("Counts preview (5 genes x 5 samples):")
print(count_mat[1:min(5,nrow(count_mat)), 1:min(5,ncol(count_mat)), drop = FALSE])

message("Library size summary (raw counts):")
print(summary(colSums(count_mat)))

message("========== [2] Filter low-expression genes (CPM rule) ==========")
libsize <- colSums(count_mat)
cpm_simple <- sweep(count_mat, 2, libsize / 1e6, "/")

keep <- rowSums(cpm_simple >= MIN_CPM) >= MIN_SAMPLES
message("Keep genes with CPM >= ", MIN_CPM, " in >= ", MIN_SAMPLES, " samples: ",
        sum(keep), " / ", nrow(count_mat))
count_mat <- count_mat[keep, , drop = FALSE]

message("After filter dim: ", paste(dim(count_mat), collapse = " x "))

message("========== [3] Compute logCPM ==========")
if (USE_EDGE_R && !requireNamespace("edgeR", quietly = TRUE)) {
  message("edgeR not installed; fallback to simple log2(CPM+1) (no TMM).")
  USE_EDGE_R <- FALSE
}

if (USE_EDGE_R) {
  y <- edgeR::DGEList(counts = count_mat)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)  # log2(CPM + 1 prior)
} else {
  libsize2 <- colSums(count_mat)
  cpm2 <- sweep(count_mat, 2, libsize2 / 1e6, "/")
  log_cpm <- log2(cpm2 + 1)
}

message("logCPM dim (genes x samples): ", paste(dim(log_cpm), collapse = " x "))
message("logCPM preview (5 genes x 5 samples):")
print(log_cpm[1:min(5,nrow(log_cpm)), 1:min(5,ncol(log_cpm)), drop = FALSE])

message("========== [4] Remove problematic genes for WGCNA ==========")
# Remove all-zero genes (should be gone already)
nonzero <- rowSums(count_mat) > 0
if (sum(!nonzero) > 0) {
  message("Dropping all-zero genes: ", sum(!nonzero))
  log_cpm <- log_cpm[nonzero, , drop = FALSE]
}

# Remove zero-variance genes (critical for correlation matrix stability)
v <- apply(log_cpm, 1, var, na.rm = TRUE)
keep_var <- v > 0
message("Zero-variance genes dropped: ", sum(!keep_var))
log_cpm <- log_cpm[keep_var, , drop = FALSE]

message("After var-filter dim: ", paste(dim(log_cpm), collapse = " x "))

# Optional: keep top variable genes
if (!is.na(TOP_VAR_GENES) && nrow(log_cpm) > TOP_VAR_GENES) {
  vv <- apply(log_cpm, 1, var)
  ord <- order(vv, decreasing = TRUE)
  log_cpm <- log_cpm[ord[1:TOP_VAR_GENES], , drop = FALSE]
  message("Keeping top variable genes: ", TOP_VAR_GENES)
}

message("Final gene count: ", nrow(log_cpm))

message("========== [5] Build WGCNA datExpr (samples x genes) ==========")
datExpr <- t(log_cpm)  # samples as rows, genes as columns
message("datExpr dim (samples x genes): ", paste(dim(datExpr), collapse = " x "))
message("datExpr preview (3 samples x 5 genes):")
print(datExpr[1:min(3,nrow(datExpr)), 1:min(5,ncol(datExpr)), drop = FALSE])

message("========== [6] Save outputs ==========")
out_datExpr <- file.path(OUT_DIR, paste0(PREFIX, "_datExpr_samplesXgenes_logCPM.csv"))
out_logcpm  <- file.path(OUT_DIR, paste0(PREFIX, "_logCPM_genesXsamples.csv"))

fwrite(as.data.table(datExpr, keep.rownames = "sample"), out_datExpr, sep = ",", quote = TRUE)
fwrite(as.data.table(log_cpm, keep.rownames = "gene_id"), out_logcpm, sep = ",", quote = TRUE)

message("Done.\nWGCNA input datExpr CSV: ", out_datExpr)

