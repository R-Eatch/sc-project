#!/usr/bin/env Rscript
# ============================================================
# WGCNA pipeline (Heart HF-LV, Lung Infection) + logging + QC exports
# Input:  datExpr CSV where rows = samples, columns = genes (numeric)
# Output: modules, module-trait, module gene lists, kME plots, top20 network exports
#
# What you asked (implemented):
#   - 默认不做高变筛选 (DO_TOP_VAR_GENES = FALSE)
#   - 全流程写入 log 文件 + 关键统计导出 (qc_summary.csv)
#   - 记录分析耗时
#   - module-trait heatmap 自动根据模块数调画布/字号，避免文字挤在一起
#   - 导出每个模块 top20 (|kME|) 的互作网络到 Cytoscape 文件（edges/nodes）
# ============================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# -----------------------------
# User settings (EDIT ME)
# -----------------------------
DATASETS <- list(
  heart = list(
    datExpr_csv = "D:/111/heartWGCNA/data/GSE161472_FPKM/preprocessed/GSE161472_ALL_preprocessed_datExpr_samplesXgenes.csv",
    outdir      = "D:/111/heartWGCNA/result/WGCNA_heart",
    trait_name  = "HF",
    pos_regex   = "(HFrEF)",
    neg_regex   = "(NF)",
    # optional: use metadata instead of regex parsing
    metadata_csv = NULL,
    trait_col    = NULL
  ),
  lung = list(
    datExpr_csv = "D:/111/heartWGCNA/data/GSE262433_RAW/step1_merge_counts/GSE262433_datExpr_samplesXgenesname_logCPM.csv",
    outdir      = "D:/111/heartWGCNA/result/WGCNA_lung",
    trait_name  = "Infection",
    pos_regex   = "(PAO1_calu3)",
    neg_regex   = "^(?!.*PAO1_calu3).*calu3",
    # optional: use metadata instead of regex parsing
    metadata_csv = NULL,
    trait_col    = NULL
  )
)

# -----------------------------
# Basic thresholds
# -----------------------------
P_THRESH   <- 0.05
COR_THRESH <- 0.30

# Soft power
POWERS <- c(1:10, seq(12, 30, by = 2))
TARGET_SFT_R2 <- 0.85

# -----------------------------
# Filtering
# -----------------------------
# 高变筛选
DO_TOP_VAR_GENES <- TRUE
TOP_VAR_N <- 12000

# 低表达过滤
DO_MEAN_FILTER <- FALSE
MEAN_LOGCPM_MIN <- 0.5

# -----------------------------
# Module detection parameters (direct tuning)
# -----------------------------
MIN_MODULE_SIZE     <- 50
DEEP_SPLIT          <- 2
MERGE_CUT_HEIGHT    <- 0.35
PAM_RESPECTS_DENDRO <- TRUE
REASSIGN_THRESHOLD  <- 0

WGCNA_PARAMS <- list(
  networkType      = "signed",
  TOMType          = "signed",
  maxBlockSize     = 20000,
  numericLabels    = TRUE,
  saveTOMs         = TRUE
)

# -----------------------------
# Network export (top20 per module)
# -----------------------------
EXPORT_NETWORKS <- TRUE
NETWORK_TOP_N        <- 20
NETWORK_TOM_THRESH   <- 0.02

# Outlier cutting (optional)
DO_OUTLIER_CUT <- FALSE
OUTLIER_CUT_HEIGHT <- NULL

# -----------------------------
# Logger + utilities
# -----------------------------
.LOG_FILE <- NULL
msg <- function(...) {
  line <- paste0("[", format(Sys.time(), "%F %T"), "] ", sprintf(...))
  cat(line, "
")
  if (!is.null(.LOG_FILE)) cat(line, "
", file = .LOG_FILE, append = TRUE)
}

safe_mkdir <- function(d) if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

init_logger <- function(outdir, tag) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  lf <- file.path(outdir, sprintf("WGCNA_%s_%s.log", tag, ts))
  .LOG_FILE <<- lf
  msg("Log file: %s", lf)
  lf
}

write_summary_kv <- function(path, kv_list) {
  dt <- data.table(
    metric = names(kv_list),
    value  = as.character(unlist(kv_list))
  )
  fwrite(dt, path)
}

write_session_info <- function(path) {
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(capture.output(sessionInfo()), con)
}

read_datExpr <- function(path) {
  dt <- fread(path)
  if (ncol(dt) < 2) stop("datExpr CSV has <2 columns: ", path)
  
  first_col <- dt[[1]]
  if (!is.numeric(first_col)) {
    rn <- as.character(first_col)
    dt[[1]] <- NULL
    mat <- as.matrix(dt)
    rownames(mat) <- rn
  } else {
    mat <- as.matrix(dt)
  }
  
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (anyNA(mat)) {
    msg("datExpr contains NA; replacing NA with 0")
    mat[is.na(mat)] <- 0
  }
  if (nrow(mat) < 2 || ncol(mat) < 2) stop("datExpr too small after load: ", path)
  mat
}

remove_bad_genes <- function(datExpr) {
  gene_var <- apply(datExpr, 2, var)
  keep <- is.finite(gene_var) & gene_var > 0
  dropped <- sum(!keep)
  if (dropped > 0) msg("Dropping %d genes with zero/NA variance", dropped)
  datExpr[, keep, drop = FALSE]
}

mean_filter <- function(datExpr, mean_min = 0.5) {
  m <- colMeans(datExpr)
  keep <- is.finite(m) & (m >= mean_min)
  dropped <- sum(!keep)
  if (dropped > 0) msg("Dropping %d genes with mean(logCPM) < %.3f", dropped, mean_min)
  datExpr[, keep, drop = FALSE]
}

top_variable_genes <- function(datExpr, n = 12000) {
  if (ncol(datExpr) <= n) return(datExpr)
  v <- apply(datExpr, 2, var)
  v[!is.finite(v)] <- -Inf
  keep <- order(v, decreasing = TRUE)[seq_len(n)]
  msg("Keeping top %d variable genes (from %d)", n, ncol(datExpr))
  datExpr[, keep, drop = FALSE]
}

make_trait_from_regex <- function(sample_ids, trait_name, pos_regex, neg_regex) {
  s <- as.character(sample_ids)
  pos <- grepl(pos_regex, s, ignore.case = TRUE, perl = TRUE)
  neg <- grepl(neg_regex, s, ignore.case = TRUE, perl = TRUE)
  
  trait <- rep(NA_integer_, length(s))
  trait[pos & !neg] <- 1L
  trait[neg & !pos] <- 0L
  
  amb <- which(pos & neg)
  if (length(amb) > 0) msg("WARNING: %d samples match BOTH pos and neg regex -> NA. Fix regex.", length(amb))
  unk <- which(!pos & !neg)
  if (length(unk) > 0) msg("WARNING: %d samples match NEITHER pos nor neg regex -> NA. Fix regex.", length(unk))
  
  df <- data.frame(row.names = s)
  df[[trait_name]] <- trait
  df
}

make_trait_from_metadata <- function(sample_ids, metadata_csv, trait_col, trait_name) {
  md <- fread(metadata_csv)
  if (!"sample_id" %in% colnames(md)) stop("metadata must contain column: sample_id")
  if (!trait_col %in% colnames(md)) stop("metadata missing trait column: ", trait_col)
  
  md$sample_id <- as.character(md$sample_id)
  s <- as.character(sample_ids)
  hit <- match(s, md$sample_id)
  trait <- md[[trait_col]][hit]
  
  if (is.character(trait) || is.factor(trait)) {
    t0 <- tolower(as.character(trait))
    trait01 <- rep(NA_integer_, length(t0))
    trait01[t0 %in% c("0", "nf", "normal", "control", "ctrl", "mock", "pbs")] <- 0L
    trait01[t0 %in% c("1", "hf", "hfre f", "hfr ef", "infection", "infected", "pao1", "pa")] <- 1L
    still <- which(is.na(trait01))
    if (length(still) > 0) {
      suppressWarnings(num <- as.numeric(t0))
      trait01[still] <- ifelse(is.finite(num[still]), as.integer(num[still]), NA_integer_)
    }
    trait <- trait01
  } else {
    suppressWarnings(trait <- as.integer(trait))
  }
  
  df <- data.frame(row.names = s)
  df[[trait_name]] <- trait
  df
}

trait_sanity_check <- function(traits, trait_name, outdir) {
  vec <- traits[[trait_name]]
  tab_all <- table(ifelse(is.na(vec), "NA", as.character(vec)))
  fwrite(data.table(level = names(tab_all), n = as.integer(tab_all)), file.path(outdir, "traits_counts.csv"))
  
  ok <- vec[!is.na(vec)]
  if (length(ok) < 4) stop("Too few samples with non-NA trait (", length(ok), "). Fix trait parsing or metadata.")
  if (length(unique(ok)) < 2) stop("Trait has only ONE class after parsing (all ", unique(ok), "). Need both groups.")
  invisible(TRUE)
}

plot_soft_threshold <- function(sft, out_png, powers) {
  png(out_png, width = 1400, height = 600, res = 150)
  par(mfrow = c(1,2))
  cex1 <- 0.9
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
       type = "n", main = "Scale independence")
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels = powers, cex = cex1, col = "red")
  abline(h = TARGET_SFT_R2, col = "red")
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = "Mean connectivity")
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
  dev.off()
}

choose_power <- function(sft, target_R2 = 0.85) {
  fi <- sft$fitIndices
  R2 <- -sign(fi[,3]) * fi[,2]
  ok <- which(R2 >= target_R2)
  if (length(ok) == 0) return(fi[which.max(R2), 1])
  fi[min(ok), 1]
}

plot_sample_clustering <- function(datExpr, out_png, cutHeight = NULL) {
  sampleTree <- hclust(dist(datExpr), method = "average")
  png(out_png, width = 1400, height = 800, res = 150)
  par(cex = 0.8)
  par(mar = c(0,4,2,0))
  plot(sampleTree, main = "Sample clustering", sub = "", xlab = "")
  if (!is.null(cutHeight)) abline(h = cutHeight, col = "red")
  dev.off()
  list(sampleTree = sampleTree)
}

module_trait_assoc <- function(MEs, traits) {
  stopifnot(all(rownames(MEs) == rownames(traits)))
  modTraitCor <- cor(MEs, traits, use = "p")
  modTraitP   <- corPvalueStudent(modTraitCor, nSamples = nrow(MEs))
  list(cor = modTraitCor, p = modTraitP)
}

plot_module_trait_heatmap <- function(modTraitCor, modTraitP, out_png) {
  # 动态调画布和字号，并减少每格文本量（cor + 星号）
  nMod <- nrow(modTraitCor)
  nTrait <- ncol(modTraitCor)
  cexTxt <- max(0.30, min(0.80, 15 / nMod))
  w <- ifelse(nTrait <= 2, 900, 1100)
  h <- max(700, 28 * nMod)
  
  star <- matrix("", nrow = nMod, ncol = nTrait)
  star[modTraitP <= 0.001] <- "***"
  star[modTraitP > 0.001 & modTraitP <= 0.01] <- "**"
  star[modTraitP > 0.01 & modTraitP <= 0.05] <- "*"
  
  textMatrix <- paste0(signif(modTraitCor, 2), star)
  dim(textMatrix) <- dim(modTraitCor)
  
  png(out_png, width = w, height = h, res = 150)
  par(mar = c(6, 10, 3, 3))
  labeledHeatmap(
    Matrix = modTraitCor,
    xLabels = colnames(modTraitCor),
    yLabels = rownames(modTraitCor),
    ySymbols = rownames(modTraitCor),
    colorLabels = FALSE,
    colors = blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = cexTxt,
    zlim = c(-1, 1),
    main = "Module-trait relationships (cor; * p<=0.05, ** p<=0.01, *** p<=0.001)"
  )
  dev.off()
}

export_module_gene_lists <- function(moduleColors, outdir) {
  dt <- data.table(gene = names(moduleColors), module = unname(moduleColors))
  fwrite(dt, file.path(outdir, "module_gene_membership_all.csv"))
  
  mods <- sort(unique(dt$module))
  safe_mkdir(file.path(outdir, "modules"))
  for (m in mods) {
    if (m == "grey") next
    genes <- dt[module == m, gene]
    fwrite(data.table(gene = genes), file.path(outdir, "modules", paste0("genes_", m, ".csv")))
  }
}

export_kME_and_plots <- function(datExpr, MEs, traits, moduleColors, outdir, trait_col) {
  trait_vec <- traits[[trait_col]]
  keep <- which(!is.na(trait_vec))
  datExpr2 <- datExpr[keep, , drop = FALSE]
  MEs2     <- MEs[keep, , drop = FALSE]
  trait2   <- trait_vec[keep]
  
  GS  <- as.numeric(cor(datExpr2, trait2, use = "p"))
  GSp <- as.numeric(corPvalueStudent(GS, nSamples = nrow(datExpr2)))
  
  kME  <- cor(datExpr2, MEs2, use = "p")
  kMEp <- corPvalueStudent(kME, nSamples = nrow(datExpr2))
  
  gene_names <- colnames(datExpr)
  geneInfo <- data.table(gene = gene_names, module = moduleColors, GS = GS, GSp = GSp)
  for (k in seq_len(ncol(kME))) {
    me <- colnames(kME)[k]
    geneInfo[[paste0("kME_", me)]]  <- kME[, k]
    geneInfo[[paste0("kMEp_", me)]] <- kMEp[, k]
  }
  fwrite(geneInfo, file.path(outdir, "geneInfo_GS_kME_all.csv"))
  
  plot_dir <- file.path(outdir, "kME_plots")
  safe_mkdir(plot_dir)
  mods <- sort(unique(moduleColors))
  for (m in mods) {
    if (m == "grey") next
    me <- paste0("ME", m)
    kcol <- paste0("kME_", me)
    if (!kcol %in% colnames(geneInfo)) next
    
    sub <- geneInfo[module == m]
    png(file.path(plot_dir, paste0("kME_hist_", m, ".png")), width = 900, height = 600, res = 120)
    hist(sub[[kcol]], main = paste0("kME distribution: ", m, " (", nrow(sub), " genes)"), xlab = kcol)
    dev.off()
    
    png(file.path(plot_dir, paste0("GS_vs_kME_", m, ".png")), width = 900, height = 700, res = 120)
    plot(sub[[kcol]], sub$GS, pch = 16, cex = 0.5,
         xlab = kcol, ylab = "GS (gene-trait cor)", main = paste0("GS vs kME: ", m))
    abline(h = 0, v = 0, lty = 2)
    dev.off()
  }
  
  geneInfo
}

export_module_networks <- function(datExpr, moduleColors, geneInfo, outdir, softPower, networkType, TOMType) {
  if (!EXPORT_NETWORKS) return(invisible(NULL))
  netdir <- file.path(outdir, "module_networks")
  safe_mkdir(netdir)
  
  mods <- sort(unique(moduleColors))
  for (m in mods) {
    if (m == "grey") next
    
    me <- paste0("ME", m)
    kcol <- paste0("kME_", me)
    if (!kcol %in% colnames(geneInfo)) next
    
    sub <- geneInfo[module == m]
    if (nrow(sub) < 10) next
    sub$absKME <- abs(sub[[kcol]])
    sub <- sub[order(-sub$absKME)]
    
    topN <- min(NETWORK_TOP_N, nrow(sub))
    sel <- sub$gene[seq_len(topN)]
    expr_sel <- datExpr[, sel, drop = FALSE]
    
    msg("Exporting network for module %s using top %d genes", m, ncol(expr_sel))
    TOM <- TOMsimilarityFromExpr(expr_sel, power = softPower, networkType = networkType, TOMType = TOMType)
    diag(TOM) <- 0
    
    mdir <- file.path(netdir, paste0("module_", m))
    safe_mkdir(mdir)
    
    nodeFile <- file.path(mdir, paste0("cyto_nodes_", m, ".txt"))
    edgeFile <- file.path(mdir, paste0("cyto_edges_", m, ".txt"))
    
    exportNetworkToCytoscape(
      TOM,
      edgeFile = edgeFile,
      nodeFile = nodeFile,
      weighted = TRUE,
      threshold = NETWORK_TOM_THRESH,
      nodeNames = colnames(expr_sel),
      altNodeNames = colnames(expr_sel),
      nodeAttr = rep(m, ncol(expr_sel))
    )
    
    if (requireNamespace("igraph", quietly = TRUE)) {
      ig <- igraph::read_graph(edgeFile, format = "ncol", directed = FALSE)
      png(file.path(mdir, paste0("igraph_", m, ".png")), width = 1100, height = 900, res = 140)
      plot(ig, vertex.size = 4, vertex.label = NA, edge.width = 1)
      dev.off()
    }
  }
  
  invisible(TRUE)
}

# -----------------------------
# Main per-dataset runner
# -----------------------------
run_wgcna_one <- function(cfg) {
  safe_mkdir(cfg$outdir)
  start_time <- Sys.time()
  init_logger(cfg$outdir, cfg$trait_name)
  
  qc <- list()
  qc$dataset <- cfg$trait_name
  qc$datExpr_csv <- cfg$datExpr_csv
  qc$start_time <- format(start_time, "%F %T")
  
  msg("=== Dataset: %s ===", cfg$trait_name)
  msg("Reading datExpr: %s", cfg$datExpr_csv)
  
  datExpr <- read_datExpr(cfg$datExpr_csv)
  qc$n_samples_raw <- nrow(datExpr)
  qc$n_genes_raw <- ncol(datExpr)
  
  qc$na_count_raw <- sum(is.na(datExpr))
  qc$nonfinite_count_raw <- sum(!is.finite(datExpr))
  
  datExpr <- remove_bad_genes(datExpr)
  qc$n_genes_after_zeroVar <- ncol(datExpr)
  
  if (DO_MEAN_FILTER) {
    n0 <- ncol(datExpr)
    datExpr <- mean_filter(datExpr, MEAN_LOGCPM_MIN)
    qc$n_genes_dropped_meanFilter <- n0 - ncol(datExpr)
  } else {
    qc$n_genes_dropped_meanFilter <- 0
  }
  qc$n_genes_after_meanFilter <- ncol(datExpr)
  
  if (DO_TOP_VAR_GENES) {
    n0 <- ncol(datExpr)
    datExpr <- top_variable_genes(datExpr, TOP_VAR_N)
    qc$n_genes_dropped_topVar <- n0 - ncol(datExpr)
  } else {
    qc$n_genes_dropped_topVar <- 0
  }
  qc$n_genes_after_topVar <- ncol(datExpr)
  
  # goodSamplesGenes
  gsg <- goodSamplesGenes(datExpr, verbose = 3)
  qc$goodSamplesGenes_allOK <- gsg$allOK
  if (!gsg$allOK) {
    msg("goodSamplesGenes flagged issues; removing offending samples/genes")
    qc$n_bad_samples <- sum(!gsg$goodSamples)
    qc$n_bad_genes <- sum(!gsg$goodGenes)
    if (sum(!gsg$goodSamples) > 0) datExpr <- datExpr[gsg$goodSamples, , drop = FALSE]
    if (sum(!gsg$goodGenes) > 0)   datExpr <- datExpr[, gsg$goodGenes, drop = FALSE]
  } else {
    qc$n_bad_samples <- 0
    qc$n_bad_genes <- 0
  }
  qc$n_samples_after_gsg <- nrow(datExpr)
  qc$n_genes_after_gsg <- ncol(datExpr)
  
  # trait
  sample_ids <- rownames(datExpr)
  if (!is.null(cfg$metadata_csv) && !is.null(cfg$trait_col)) {
    traits <- make_trait_from_metadata(sample_ids, cfg$metadata_csv, cfg$trait_col, cfg$trait_name)
  } else {
    traits <- make_trait_from_regex(sample_ids, cfg$trait_name, cfg$pos_regex, cfg$neg_regex)
  }
  
  fwrite(data.table(sample_id = rownames(traits), trait = traits[[cfg$trait_name]]),
         file.path(cfg$outdir, "traits_parsed.csv"))
  
  trait_sanity_check(traits, cfg$trait_name, cfg$outdir)
  tv <- traits[[cfg$trait_name]]
  qc$trait_nonNA_n <- sum(!is.na(tv))
  qc$trait_NA_n <- sum(is.na(tv))
  qc$trait_0_n <- sum(tv == 0, na.rm = TRUE)
  qc$trait_1_n <- sum(tv == 1, na.rm = TRUE)
  
  # sample clustering
  sampleTree <- plot_sample_clustering(datExpr, file.path(cfg$outdir, "01_sample_clustering.png"))$sampleTree
  if (DO_OUTLIER_CUT) {
    cutH <- OUTLIER_CUT_HEIGHT
    if (is.null(cutH)) cutH <- 0.9 * max(sampleTree$height)
    clust <- cutreeStatic(sampleTree, cutHeight = cutH, minSize = 2)
    keep <- (clust == 1)
    msg("Outlier cut enabled: keeping %d/%d samples", sum(keep), length(keep))
    datExpr <- datExpr[keep, , drop = FALSE]
    traits  <- traits[rownames(datExpr), , drop = FALSE]
    plot_sample_clustering(datExpr, file.path(cfg$outdir, "01b_sample_clustering_after_cut.png"), cutHeight = cutH)
    qc$outlier_cut_height <- cutH
  } else {
    qc$outlier_cut_height <- NA
  }
  qc$n_samples_final <- nrow(datExpr)
  qc$n_genes_final <- ncol(datExpr)
  
  # soft threshold
  msg("Picking soft threshold...")
  sft <- pickSoftThreshold(datExpr, powerVector = POWERS, networkType = WGCNA_PARAMS$networkType, verbose = 5)
  plot_soft_threshold(sft, file.path(cfg$outdir, "02_soft_threshold.png"), POWERS)
  softPower <- choose_power(sft, target_R2 = TARGET_SFT_R2)
  qc$softPower <- softPower
  msg("Selected soft power = %d", softPower)
  
  # modules
  msg("Running blockwiseModules... (minModuleSize=%d, deepSplit=%d, mergeCutHeight=%.2f)",
      MIN_MODULE_SIZE, DEEP_SPLIT, MERGE_CUT_HEIGHT)
  tomBase <- file.path(cfg$outdir, "TOM")
  
  net <- blockwiseModules(
    datExpr,
    power = softPower,
    networkType = WGCNA_PARAMS$networkType,
    TOMType = WGCNA_PARAMS$TOMType,
    minModuleSize = MIN_MODULE_SIZE,
    reassignThreshold = REASSIGN_THRESHOLD,
    mergeCutHeight = MERGE_CUT_HEIGHT,
    deepSplit = DEEP_SPLIT,
    pamRespectsDendro = PAM_RESPECTS_DENDRO,
    maxBlockSize = WGCNA_PARAMS$maxBlockSize,
    numericLabels = WGCNA_PARAMS$numericLabels,
    saveTOMs = WGCNA_PARAMS$saveTOMs,
    saveTOMFileBase = tomBase,
    verbose = 5
  )
  
  moduleColors <- labels2colors(net$colors)
  names(moduleColors) <- colnames(datExpr)
  MEs <- orderMEs(net$MEs)
  labels <- net$colors
  cols   <- labels2colors(labels)
  
  map <- data.table(
    module_label = as.integer(sort(unique(labels))),
    module_color = labels2colors(sort(unique(labels)))
  )
  
  # 模块大小（按 label）
  size_label <- as.data.table(table(labels))
  setnames(size_label, c("labels","N"), c("module_label","n_genes"))
  size_label[, module_label := as.integer(as.character(module_label))]
  map <- merge(map, size_label, by="module_label", all.x=TRUE)
  fwrite(map[order(module_label)], file.path(cfg$outdir, "module_label_color_map.csv"))
  
  # dendrograms
  safe_mkdir(file.path(cfg$outdir, "dendrograms"))
  for (b in seq_along(net$dendrograms)) {
    png(file.path(cfg$outdir, "dendrograms", sprintf("03_dendrogram_block%02d.png", b)),
        width = 1600, height = 850, res = 150)
    plotDendroAndColors(
      net$dendrograms[[b]],
      labels2colors(net$colors[net$blockGenes[[b]]]),
      "Module colors",
      dendroLabels = FALSE,
      hang = 0.03,
      addGuide = TRUE,
      guideHang = 0.05,
      main = paste0("Gene dendrogram - block ", b)
    )
    dev.off()
  }
  
  # module-trait
  trait_vec <- traits[[cfg$trait_name]]
  keep_trait <- which(!is.na(trait_vec))
  assoc <- module_trait_assoc(MEs[keep_trait, , drop = FALSE], traits[keep_trait, , drop = FALSE])
  fwrite(as.data.table(assoc$cor, keep.rownames = "module"), file.path(cfg$outdir, "04_module_trait_cor.csv"))
  fwrite(as.data.table(assoc$p,   keep.rownames = "module"), file.path(cfg$outdir, "04_module_trait_p.csv"))
  plot_module_trait_heatmap(assoc$cor, assoc$p, file.path(cfg$outdir, "04_module_trait_heatmap.png"))
  
  cor_vec <- assoc$cor[, cfg$trait_name]
  p_vec   <- assoc$p[, cfg$trait_name]
  sig <- which(is.finite(cor_vec) & is.finite(p_vec) & (abs(cor_vec) >= COR_THRESH) & (p_vec <= P_THRESH))
  sig_modules <- rownames(assoc$cor)[sig]
  fwrite(data.table(module = sig_modules, cor = cor_vec[sig], p = p_vec[sig]),
         file.path(cfg$outdir, "05_significant_modules.csv"))
  
  qc$n_modules_total_exclGrey <- length(setdiff(unique(moduleColors), "grey"))
  qc$n_sig_modules <- length(sig_modules)
  
  # exports
  export_module_gene_lists(moduleColors, cfg$outdir)
  geneInfo <- export_kME_and_plots(datExpr, MEs, traits, moduleColors, cfg$outdir, cfg$trait_name)
  export_module_networks(datExpr, moduleColors, geneInfo, cfg$outdir, softPower, WGCNA_PARAMS$networkType, WGCNA_PARAMS$TOMType)
  
  # module size stats
  mod_sizes <- sort(table(moduleColors), decreasing = TRUE)
  fwrite(data.table(module = names(mod_sizes), n_genes = as.integer(mod_sizes)),
         file.path(cfg$outdir, "module_sizes.csv"))
  
  # save rds
  saveRDS(list(
    datExpr_path  = cfg$datExpr_csv,
    datExpr_dim   = dim(datExpr),
    traits        = traits,
    softPower     = softPower,
    params        = list(
      MIN_MODULE_SIZE = MIN_MODULE_SIZE,
      DEEP_SPLIT = DEEP_SPLIT,
      MERGE_CUT_HEIGHT = MERGE_CUT_HEIGHT,
      PAM_RESPECTS_DENDRO = PAM_RESPECTS_DENDRO,
      REASSIGN_THRESHOLD = REASSIGN_THRESHOLD,
      networkType = WGCNA_PARAMS$networkType,
      TOMType = WGCNA_PARAMS$TOMType,
      DO_TOP_VAR_GENES = DO_TOP_VAR_GENES,
      TOP_VAR_N = TOP_VAR_N,
      DO_MEAN_FILTER = DO_MEAN_FILTER,
      MEAN_LOGCPM_MIN = MEAN_LOGCPM_MIN
    ),
    moduleColors  = moduleColors,
    MEs           = MEs,
    moduleTraitCor= assoc$cor,
    moduleTraitP  = assoc$p,
    significantModules = sig_modules,
    geneInfo_path = file.path(cfg$outdir, "geneInfo_GS_kME_all.csv")
  ), file = file.path(cfg$outdir, "WGCNA_result.rds"))
  
  # qc + runtime
  end_time <- Sys.time()
  qc$end_time <- format(end_time, "%F %T")
  qc$elapsed_minutes <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 3)
  write_summary_kv(file.path(cfg$outdir, "qc_summary.csv"), qc)
  write_session_info(file.path(cfg$outdir, "sessionInfo.txt"))
  
  msg("QC summary saved: %s", file.path(cfg$outdir, "qc_summary.csv"))
  msg("Done. Outputs in: %s", cfg$outdir)
  msg("Total elapsed: %.2f minutes", qc$elapsed_minutes)
  
  invisible(TRUE)
}

# -----------------------------
# Run all datasets
# -----------------------------
msg("Starting WGCNA for %d datasets", length(DATASETS))
for (nm in names(DATASETS)) {
  run_wgcna_one(DATASETS[[nm]])
}
msg("All done.")
