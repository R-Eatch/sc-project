#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# hdWGCNA_v5_unsplit_symbol_pipeline.R
# Goal: work from an UNSPLIT Seurat object; create age+celltype labels; subset by a chosen age-celltype;
#       build an RNA_sym assay with SYMBOL rownames, WITHOUT re-running Normalize/UMAP/PCA, etc.
#       Reuse existing normalized data if available; otherwise, stop with a clear message.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(Matrix)
  library(clusterProfiler)
  library(org.Mm.eg.db)  # mouse
})
suppressPackageStartupMessages({ library(future) })

options(future.globals.maxSize = 200 * 1024^3)

# 避免 Linux 下 fork 复制巨型对象
Sys.setenv(R_FUTURE_FORK_ENABLE = "FALSE")

enableWGCNAThreads(nThreads = 4)
theme_set(theme_cowplot())
set.seed(12345)

# ===== User-configurable variables =====
# Paths
INPUT_RDS       <- "../../../data/splits_seurat/combined_seurat_full.rds"   # <-- unsplit Seurat object
result_root     <- "./result_unsplit_hdWGCNA"                 # root output dir

# Metadata columns present in the UNSPLIT object
AGE_COL         <- "Age_group"               # e.g., "03_months", "06_months", ...
TYPE_COL        <- "Main_cell_type"    # e.g., "Ventricular cardiomyocytes"
GROUPBY         <- "age_plus_celltype"            # for metacells & summaries

# Choose ONE age + ONE cell type for analysis in this run
TARGET_AGE         <- "23_months"
TARGET_CELLTYPE    <- "Fibroblasts"

# Dataset tags for filenames/titles
DATASET            <- "H"
USE_RANDOM_SUBSET  <- TRUE
RUN_KEGG <- FALSE
AGE_TAG            <- ""
SPEFIC_CELLTYPE    <- paste(TARGET_AGE, TARGET_CELLTYPE, sep = "+")

# hdWGCNA parameters (kept same as your prior script where possible)
soft_power      <- NULL
deepSplit       <- 4
detectCutHeight <- 0.995
minModuleSize   <- 10
mergeCutHeight  <- 0.1
kcells          <- 30
n_genes         <- 5
SUBSET_NUMBER   <- 50000
# ===== helpers =====
PREFIX <- paste(DATASET, AGE_TAG, SPEFIC_CELLTYPE, sep = "-")
msg_ok <- function(...) message("✅ ", sprintf(...))
msg_sk <- function(...) message("⏭️ ", sprintf(...))
msg_er <- function(...) stop(sprintf(...))

plot_module_native <- function(obj, resultpath, dataset, age_tag, spefic_celltype,
                               n_genes=5, group_var=GROUPBY,
                               reduction="umap", assay="RNA_sym") {
  modules_tbl <- GetModules(obj) %>%
    dplyr::filter(module != "grey") %>%
    dplyr::select(gene_name, module, dplyr::starts_with("kME_"))
  if (nrow(modules_tbl) == 0) return(invisible(NULL))
  hub_tbl <- modules_tbl %>%
    dplyr::group_by(module) %>%
    dplyr::group_modify(~{
      col_kme <- paste0("kME_", .y$module)
      .x %>% dplyr::arrange(dplyr::desc(.data[[col_kme]])) %>% dplyr::slice_head(n = n_genes)
    }) %>% dplyr::ungroup()

  mod_list <- split(hub_tbl, hub_tbl$module)
  fp_list <- list()
  for (mod in names(mod_list)) {
    genes_mod <- mod_list[[mod]]$gene_name
    genes_mod <- intersect(genes_mod, rownames(obj[[assay]]))
    if (!length(genes_mod)) next
    fp_list[[mod]] <- FeaturePlot(obj, features=genes_mod, reduction=reduction,
                                  ncol=length(genes_mod), order=TRUE) +
      plot_annotation(title=paste0(mod, " hub genes: ", paste(genes_mod, collapse=", "))) +
      theme(plot.title = element_text(hjust=.5))
  }
  if (length(fp_list)) {
    p_feature <- patchwork::wrap_plots(fp_list, ncol=1)
    ggsave(file.path(resultpath, paste0(PREFIX, "-FeaturePlot.png")),
           p_feature, width=15, height=2 + 2.5*length(fp_list), dpi=300, bg='white', limitsize=FALSE)
  }
  gene_order <- unique(hub_tbl$gene_name)
  gene_order <- intersect(gene_order, rownames(obj[[assay]]))
  if (length(gene_order)) {
    p_dot <- DotPlot(obj, features=gene_order, group.by=group_var, assay=assay) +
      labs(y="Gene (module)", x=group_var,
           title=paste0("DotPlot of hub genes — ", dataset, " / ", age_tag, " / ", spefic_celltype)) +
      theme(axis.text.y = element_text(size=8))
    ggsave(file.path(resultpath, paste0(PREFIX, "-DotPlot.png")),
           p_dot, width=10 + 2.5*length(unique(hub_tbl$module)), height=6, dpi=300, bg='white', limitsize=FALSE)
  }
  invisible(TRUE)
}

module_enrich <- function(modules_df, result_path, dataset, age_tag, spefic_celltype,
                          gene_col="gene_name", cor_threshold=0.5, min_gene_count=3,
                          top_show=15, plot_width=14, plot_height=14, RUN_KEGG=TRUE) {

  OrgDb    <- org.Mm.eg.db
  kegg_org <- "mmu"

  plot_dir  <- file.path(result_path, "plots");  dir.create(plot_dir,  TRUE, TRUE)
  table_dir <- file.path(result_path, "tables"); dir.create(table_dir, TRUE, TRUE)

  # 占位空白图（完全留白）
  blank_plot <- function() ggplot() + theme_void()

  # 占位提示图（给 GO 用，便于查错）
  empty_plot <- function(txt) ggplot() + annotate("text", x=.5, y=.5, label=txt, size=4) + theme_void()

  save_tbl   <- function(obj, path) write.csv(as.data.frame(obj), path, row.names = FALSE)

  # 安全 dotplot：判空、限制 showCategory、清 NA，并在构建期加 tryCatch
  safe_dotplot <- function(x, top_show = 15, main_title = "", blank_when_empty = FALSE) {
    if (is.null(x)) return(if (blank_when_empty) blank_plot() else empty_plot(paste(main_title, "no result")))
    df <- tryCatch(as.data.frame(x), error=function(e) NULL)
    if (is.null(df) || nrow(df) == 0) return(if (blank_when_empty) blank_plot() else empty_plot(paste(main_title, "no result")))
    # 清掉常见 NA 行
    df <- df[stats::complete.cases(df), , drop = FALSE]
    if ("p.adjust" %in% names(df)) df <- df[!is.na(df$p.adjust), , drop = FALSE]
    if (nrow(df) == 0) return(if (blank_when_empty) blank_plot() else empty_plot(paste(main_title, "no significant terms")))
    show_n <- max(1, min(top_show, nrow(df)))

    # 重新构建一个 enrichResult，防止 as.data.frame 后丢失 S4（enrichplot 接受 enrichResult 或 data.frame 均可，但更稳妥）
    p <- tryCatch({
      p0 <- enrichplot::dotplot(x, showCategory = show_n) +
            ggtitle(main_title) + theme_classic(base_size = 11)
      # 关键：强制构建，提前捕获惰性期错误
      ggplot2::ggplot_build(p0)
      p0
    }, error = function(e) {
      if (blank_when_empty) {
        blank_plot()
      } else {
        empty_plot(paste(main_title, "plot error:\n", conditionMessage(e)))
      }
    })
    return(p)
  }

  all_modules <- unique(modules_df$module)

  for (mod in all_modules) {
    kme_col <- paste0("kME_", mod)
    this_df <- modules_df %>% dplyr::filter(module == !!mod)
    if (!kme_col %in% names(this_df)) next

    genes <- this_df %>%
      dplyr::filter(abs(.data[[kme_col]]) >= cor_threshold) %>%
      dplyr::pull(.data[[gene_col]]) %>% unique() %>% stats::na.omit()

    if (length(genes) < min_gene_count) next

    # ---------- GO ----------
    go_plots <- vector("list", 3); names(go_plots) <- c("BP","CC","MF")
    for (ont in c("BP","CC","MF")) {
      go_res <- tryCatch(
        enrichGO(gene = genes, OrgDb = OrgDb, keyType = "SYMBOL", ont = ont,
                 pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20),
        error = function(e) NULL
      )
      if (!is.null(go_res) && nrow(go_res@result) > 0) {
        save_tbl(go_res@result, file.path(
          table_dir, sprintf("%s-%s-%s_GO_%s.csv", dataset, age_tag, gsub('[/ ]','_', spefic_celltype), ont)
        ))
      }
      go_plots[[ont]] <- safe_dotplot(
        go_res, top_show = top_show,
        main_title = sprintf("%s / %s / %s — GO-%s", dataset, age_tag, mod, ont),
        blank_when_empty = FALSE  # GO 为空时给提示文字，便于诊断
      )
    }

    # ---------- KEGG ----------
    kegg_plot <- NULL
    if (!isTRUE(RUN_KEGG)) {
      # 需求：关闭时右下角留白
      kegg_plot <- blank_plot()
    } else {
      entrez_map <- suppressMessages(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb))
      if (is.null(entrez_map) || nrow(entrez_map) == 0) {
        kegg_plot <- empty_plot(paste(mod, "no valid genes for KEGG"))
      } else {
        kegg_res <- tryCatch(
          enrichKEGG(gene = unique(entrez_map$ENTREZID), organism = kegg_org, keyType = "ncbi-geneid",
                     pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.25),
          error = function(e) NULL
        )
        kegg_df <- tryCatch(if (!is.null(kegg_res)) as.data.frame(kegg_res) else NULL, error=function(e) NULL)
        if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
          # 把 geneID 从 ENTREZ 映射回 SYMBOL（可读性更好）
          kegg_df$geneID <- vapply(kegg_df$geneID, function(x) {
            ids  <- strsplit(x, "/")[[1]]
            syms <- entrez_map$SYMBOL[match(ids, entrez_map$ENTREZID)]
            paste(stats::na.omit(syms), collapse="/")
          }, character(1))
          save_tbl(kegg_df, file.path(
            table_dir, sprintf("%s-%s-%s_KEGG.csv", dataset, age_tag, gsub('[/ ]','_', spefic_celltype))
          ))
          kegg_plot <- safe_dotplot(
            kegg_res, top_show = top_show,
            main_title = sprintf("%s / %s / %s — KEGG", dataset, age_tag, mod),
            blank_when_empty = FALSE
          )
        } else {
          kegg_plot <- empty_plot(paste(mod, "no significant KEGG"))
        }
      }
    }

    # ---------- 合图 & 保存 ----------
    combined <- (go_plots[["BP"]] | go_plots[["CC"]]) / (go_plots[["MF"]] | kegg_plot)

    # 保存前再“试构建”一次，避免 ggsave 阶段才报错
    ok_to_save <- TRUE
    tryCatch({ ggplot2::ggplot_build(combined) }, error = function(e) { ok_to_save <<- FALSE })

    pdf_path <- file.path(plot_dir, sprintf("%s-%s-%s_GO_KEGG_%s.pdf",
                                            dataset, age_tag, gsub('[/ ]','_', spefic_celltype), mod))
    if (ok_to_save) {
      tryCatch(
        ggsave(pdf_path, combined, width = plot_width, height = plot_height, device = cairo_pdf),
        error = function(e) {
          # 兜底：cairo 不可用或保存失败则转 PNG
          ggsave(fs::path_ext_set(pdf_path, "png"), combined,
                 width = plot_width, height = plot_height, dpi = 200, bg = "white")
        }
      )
    } else {
      # 构建失败也给个兜底 PNG（用 cowplot 把错误消息写进去）
      fallback <- empty_plot("Plot build failed; saved placeholder instead")
      ggsave(fs::path_ext_set(pdf_path, "png"), fallback,
             width = plot_width, height = plot_height, dpi = 200, bg = "white")
    }
  }

  msg_ok("Enrichment finished → %s", result_path)
  invisible(TRUE)
}


# ===== Load UNSPLIT RDS & create composite age+celltype labels =====
seurat_all <- readRDS(INPUT_RDS)
seurat_all <- SeuratObject::UpdateSeuratObject(seurat_all)
if(USE_RANDOM_SUBSET){
  set.seed(123)  # 保证可复现
  cells_use <- sample(colnames(seurat_all), SUBSET_NUMBER)
  seurat_obj <- subset(seurat_all, cells = cells_use)
  seurat_all <- seurat_obj
  rm(seurat_obj)

}


meta <- seurat_all@meta.data
#if (!all(c(AGE_COL, TYPE_COL) %in% colnames(meta))) {
#  msg_er("Missing required metadata columns: need %s ; have %s",
#         paste(c(AGE_COL, TYPE_COL), collapse=", "), paste(colnames(meta), collapse=", "))
#}

# Create two convenient composite labels
meta[["age_celltype"]]      <- paste(meta[[AGE_COL]], meta[[TYPE_COL]], sep = "-")  # age-celltype
meta[["age_plus_celltype"]] <- paste(meta[[AGE_COL]], meta[[TYPE_COL]], sep = "+")  # age+celltype
seurat_all@meta.data <- meta

# Subset to the selected age & cell type
#cells_keep <- rownames(meta)[meta[[AGE_COL]] == TARGET_AGE & meta[[TYPE_COL]] == TARGET_CELLTYPE]
#if (length(cells_keep) == 0) msg_er("No cells found for %s + %s", TARGET_AGE, TARGET_CELLTYPE)
seurat_obj <- seurat_all
rm(seurat_all); invisible(gc())

# ===== Build v5 RNA_sym assay with SYMBOL rownames, reusing existing normalized data (no re-Normalize) =====
DefaultAssay(seurat_obj) <- "RNA"
cts <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
mf  <- seurat_obj[["RNA"]]@meta.features
if (!"gene_name" %in% colnames(mf)) msg_er("meta.features for RNA lacks 'gene_name'.")

# derive SYMBOLs in the SAME order as counts rows
sym <- mf[rownames(cts), "gene_name", drop=TRUE]
sym <- as.character(sym)
sym <- trimws(sym)
fallback_ids <- rownames(cts)
mask <- is.na(sym) | sym == "" | sym == "NA"
sym[mask] <- fallback_ids[mask]
sym <- trimws(sym)
mask2 <- is.na(sym) | sym == ""
sym[mask2] <- fallback_ids[mask2]
stopifnot(length(sym) == nrow(cts))
sym_unique <- make.unique(sym)

# also remap the normalized DATA slot if present — we DO NOT recompute NormalizeData
rna_data <- tryCatch(GetAssayData(seurat_obj, assay="RNA", slot="data"), error=function(e) NULL)
if (is.null(rna_data) || nrow(rna_data) != nrow(cts)) {
  msg_er("No usable normalized 'data' slot found to reuse. Please normalize upstream, or provide an assay with a populated 'data' slot.")
}
sym_clean <- gsub("_", "-", sym_unique, fixed = TRUE)
# apply SYMBOL rownames consistently
rownames(cts)      <- sym_clean
rownames(rna_data) <- sym_clean
options(Seurat.object.assay.version = "v5")
seurat_obj[["RNA_sym"]] <- CreateAssayObject(counts = cts)
# inject the reused normalized data into RNA_sym
seurat_obj[["RNA_sym"]]@data <- rna_data[rownames(seurat_obj[["RNA_sym"]]@counts), colnames(seurat_obj)]
DefaultAssay(seurat_obj) <- "RNA_sym"

# Optional quick overview (reuses existing UMAP/PCA if already present; we DO NOT recompute)
if ("umap" %in% names(seurat_obj@reductions)) {
  p_over <- DimPlot(seurat_obj, reduction="umap", group.by=GROUPBY, label=TRUE) +
    ggtitle(paste0("UMAP — ", DATASET, " / ", AGE_TAG, " / ", SPEFIC_CELLTYPE)) + NoLegend()
} else if ("pca" %in% names(seurat_obj@reductions)) {
  p_over <- DimPlot(seurat_obj, reduction="pca", group.by=GROUPBY, label=TRUE) +
    ggtitle(paste0("PCA — ", DATASET, " / ", AGE_TAG, " / ", SPEFIC_CELLTYPE)) + NoLegend()
} else {
  p_over <- ggplot() + annotate("text", x=.5, y=.5, label="No UMAP/PCA present — skip overview", size=5) + theme_void()
}
resultpath <- file.path(result_root, sprintf("result_%s_%s_%s", DATASET, gsub("[^0-9A-Za-z]+","", AGE_TAG), gsub("[^0-9A-Za-z]+","", SPEFIC_CELLTYPE)))
dir.create(resultpath, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(resultpath, paste0(PREFIX, "-overview.png")), p_over, width=7, height=6, dpi=300, bg="white")

# ===== hdWGCNA (no re-Normalize/UMAP/PCA). We will reuse existing reductions ONLY. =====
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select="fraction", fraction=0.10, wgcna_name="tutorial")

# Choose a reduction already present for metacells; prefer PCA>UMAP>tSNE>LSI
redcandidates <- c("umap","pca","tsne","lsi")
red_use <- intersect(redcandidates, names(seurat_obj@reductions))
if (length(red_use) == 0) msg_er("No existing reductions (PCA/UMAP/tSNE/LSI) found. Please compute one upstream.")
red_use <- 'pca'
msg_ok("Using existing reduction for metacells: %s", red_use)

seurat_obj <- MetacellsByGroups(seurat_obj, group.by=GROUPBY, reduction=red_use, k=kcells, max_shared=10, ident.group=GROUPBY)
seurat_obj <- NormalizeMetacells(seurat_obj)  # normalization over metacells (kept; does not touch raw cells)
seurat_obj <- SetDatExpr(seurat_obj, group_name=SPEFIC_CELLTYPE, group.by=GROUPBY, assay="RNA_sym", layer="data")

seurat_obj <- TestSoftPowers(seurat_obj, networkType="signed")
plot_list <- PlotSoftPowers(seurat_obj); wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(resultpath, paste0(PREFIX, "-softpower_table.csv")), row.names=FALSE)

seurat_obj <- ConstructNetwork(seurat_obj,
  tom_name=SPEFIC_CELLTYPE, overwrite_tom=TRUE, soft_power=soft_power,
  deepSplit=deepSplit, detectCutHeight=detectCutHeight,
  minModuleSize=minModuleSize, mergeCutHeight=mergeCutHeight)

png(file.path(resultpath, paste0(PREFIX, "-dendrogram.png")), width=2400, height=1800, res=300)
PlotDendrogram(seurat_obj, main=paste0(PREFIX, " hdWGCNA Dendrogram"))
dev.off()

TOM <- GetTOM(seurat_obj)

# No re-Scale/FindVariableFeatures/RunPCA here on single cells; only module ops
seurat_obj <- ModuleEigengenes(seurat_obj)
hMEs <- GetMEs(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj, group.by=GROUPBY, group_name=SPEFIC_CELLTYPE)

seurat_obj <- ResetModuleNames(seurat_obj, new_name=paste0(DATASET, "-", AGE_TAG, "-", SPEFIC_CELLTYPE, "-M"))
saveRDS(seurat_obj, file=file.path(resultpath, paste0(PREFIX, "-hdWGCNA_object.rds")))

modules <- GetModules(seurat_obj) %>% dplyr::filter(module != "grey")
if (nrow(modules) > 0) {
  p2 <- PlotKMEs(seurat_obj, ncol=5, n_hubs=15)
  ggsave(filename=file.path(resultpath, paste0(PREFIX, "-KMEsplot.png")),
         plot=p2, dpi=300, bg='white', width=30, height=15)
  write.csv(modules, file=file.path(resultpath, paste0(PREFIX, "-modules.csv")), row.names=FALSE)

  ModuleNetworkPlot(seurat_obj, outdir=file.path(resultpath, "moduleNetwork-default"))
  ModuleNetworkPlot(seurat_obj, outdir=file.path(resultpath, "moduleNetwork-20"),
                    n_inner=20, n_outer=30, n_conns=Inf, plot_size=c(10,10), vertex.label.cex=1)

  plot_list <- ModuleFeaturePlot(seurat_obj, features="hMEs", order=TRUE)
  p3 <- wrap_plots(plot_list, ncol=5)
  ggsave(filename=file.path(resultpath, paste0(PREFIX, "-hMEs-fp.png")),
         plot=p3, dpi=200, bg='white', width=30, height=max(6, length(unique(modules$module))*1), limitsize=FALSE)

  mods <- setdiff(levels(modules$module), "grey")
  if (length(mods)) {
    png(file.path(resultpath, paste0(PREFIX, "-HubGeneNetwork.png")), width=2400, height=1800, res=300)
    HubGeneNetworkPlot(seurat_obj, n_hubs=10, n_other=20, edge_prop=0.75, mods=mods)
    dev.off()
  }

  # Module UMAP uses module graph embedding; does not require cell UMAP recomputation
  seurat_obj <- RunModuleUMAP(seurat_obj, n_hubs=10, n_neighbors=15, min_dist=0.1)
  png(file.path(resultpath, paste0(PREFIX, "-module-umap.png")), width=1800, height=1800, res=300)
  ModuleUMAPPlot(seurat_obj, edge.alpha=0.25, sample_edges=TRUE, edge_prop=0.1, label_hubs=2, keep_grey_edges=FALSE)
  dev.off()

  plot_module_native(obj=seurat_obj, resultpath=resultpath, dataset=DATASET, age_tag=AGE_TAG,
                     spefic_celltype=SPEFIC_CELLTYPE, group_var=GROUPBY, n_genes=n_genes, assay="RNA_sym")

  module_enrich(modules, file.path(resultpath, "enrich"),
                dataset=DATASET, age_tag=AGE_TAG, spefic_celltype=SPEFIC_CELLTYPE,
                cor_threshold=0.3, top_show=15, RUN_KEGG=RUN_KEGG)
} else {
  msg_sk("No non-grey modules detected; skip plotting/enrichment.")
}

msg_ok("All done → %s", resultpath)

