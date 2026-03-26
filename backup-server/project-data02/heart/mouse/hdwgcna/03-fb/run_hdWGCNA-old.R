#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# hdWGCNA_v5_symbol_pipeline.R  (fixed: coerce gene_name to character before fallback)

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

options(future.globals.maxSize = 10 * 1024^3)
enableWGCNAThreads(nThreads = 4)
theme_set(theme_cowplot())
set.seed(12345)

# ===== Global variables (per-group) =====
dataset         <- "Heart"
AGE_TAG         <- "03_months"
spefic_celltype <- "Fibroblasts"
groupby         <- "Main_cell_type"
RDSpath         <- "../../../data/splits_seurat/seurat_03_months_fibroblasts.rds"
resultpath      <- "./result_Heart_03m_Fibro"

dir.create(resultpath, recursive = TRUE, showWarnings = FALSE)

# hdWGCNA parameters
soft_power     <- NULL
deepSplit       <- 4
detectCutHeight <- 0.98
minModuleSize   <- 15
mergeCutHeight  <- 0.3
kcells          <- 30
n_genes         <- 5

# ===== helpers =====
PREFIX <- paste(dataset, AGE_TAG, spefic_celltype, sep = "-")
msg_ok <- function(...) message("✅ ", sprintf(...))
msg_sk <- function(...) message("⏭️ ", sprintf(...))

plot_module_native <- function(obj, resultpath, dataset, age_tag, spefic_celltype,
                               n_genes=5, group_var="Main_cell_type",
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
                          top_show=15, plot_width=14, plot_height=14) {
  OrgDb    <- org.Mm.eg.db
  kegg_org <- "mmu"
  plot_dir  <- file.path(result_path, "plots");  dir.create(plot_dir,  TRUE, TRUE)
  table_dir <- file.path(result_path, "tables"); dir.create(table_dir, TRUE, TRUE)
  empty_plot <- function(txt) ggplot() + annotate("text", x=.5, y=.5, label=txt, size=4) + theme_void()
  save_tbl   <- function(obj, path) write.csv(as.data.frame(obj), path, row.names = FALSE)
  all_modules <- unique(modules_df$module)

  for (mod in all_modules) {
    kme_col <- paste0("kME_", mod)
    this_df <- modules_df %>% dplyr::filter(module == !!mod)
    if (!kme_col %in% names(this_df)) next
    genes <- this_df %>%
      dplyr::filter(abs(.data[[kme_col]]) >= cor_threshold) %>%
      dplyr::pull(.data[[gene_col]]) %>% unique() %>% na.omit()
    if (length(genes) < min_gene_count) next

    go_plots <- vector("list", 3); names(go_plots) <- c("BP","CC","MF")
    for (ont in c("BP","CC","MF")) {
      go_res <- tryCatch(
        enrichGO(gene=genes, OrgDb=OrgDb, keyType="SYMBOL", ont=ont,
                 pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.20),
        error=function(e) NULL
      )
      if (!is.null(go_res) && nrow(go_res@result) > 0) {
        save_tbl(go_res@result, file.path(table_dir,
                  sprintf("%s-%s-%s_GO_%s.csv", dataset, age_tag, gsub('[/ ]','_',spefic_celltype), ont)))
        go_plots[[ont]] <- tryCatch(
          dotplot(go_res, showCategory=top_show) +
            ggtitle(sprintf("%s / %s / %s — GO-%s", dataset, age_tag, mod, ont)) +
            theme_classic(base_size=11),
          error=function(e) empty_plot(paste(mod, "GO-", ont, "plot error"))
        )
      } else {
        go_plots[[ont]] <- empty_plot(paste(mod, "no significant\nGO", ont))
      }
    }

    entrez_map <- suppressMessages(bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb))
    if (is.null(entrez_map) || nrow(entrez_map) == 0) {
      kegg_plot <- empty_plot(paste(mod, "no valid genes for KEGG"))
    } else {
      kegg_res <- tryCatch(
        enrichKEGG(gene=unique(entrez_map$ENTREZID), organism=kegg_org, keyType="ncbi-geneid",
                   pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.25),
        error=function(e) NULL
      )
      kegg_df <- tryCatch(if (!is.null(kegg_res)) as.data.frame(kegg_res) else NULL, error=function(e) NULL)
      if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
        kegg_df$geneID <- sapply(kegg_df$geneID, function(x) {
          ids  <- strsplit(x, "/")[[1]]
          syms <- entrez_map$SYMBOL[match(ids, entrez_map$ENTREZID)]
          paste(na.omit(syms), collapse="/")
        })
        save_tbl(kegg_df, file.path(table_dir,
                  sprintf("%s-%s-%s_KEGG.csv", dataset, age_tag, gsub('[/ ]','_',spefic_celltype))))
        kegg_plot <- tryCatch(
          dotplot(kegg_res, showCategory=top_show) +
            ggtitle(sprintf("%s / %s / %s — KEGG", dataset, age_tag, mod)) +
            theme_classic(base_size=11),
          error=function(e) empty_plot(paste(mod, "KEGG plot error"))
        )
      } else {
        kegg_plot <- empty_plot(paste(mod, "no significant\nKEGG"))
      }
    }

    combined <- (go_plots[["BP"]] | go_plots[["CC"]]) / (go_plots[["MF"]] | kegg_plot)
    pdf_path <- file.path(plot_dir, sprintf("%s-%s-%s_GO_KEGG_%s.pdf",
                                            dataset, age_tag, gsub('[/ ]','_',spefic_celltype), mod))
    ggsave(pdf_path, combined, width=14, height=14, device=cairo_pdf)
  }
  msg_ok("Enrichment finished → %s", result_path)
  invisible(TRUE)
}

# ===== Load RDS & build v5 assay with symbol rownames =====
seurat_obj <- readRDS(RDSpath)
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)

DefaultAssay(seurat_obj) <- "RNA"
cts <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
mf  <- seurat_obj[["RNA"]]@meta.features
if (!"gene_name" %in% colnames(mf)) stop("meta.features for RNA lacks 'gene_name'.")

# -- FIXED: coerce to character *before* fallback & make.unique --
sym <- mf[rownames(cts), "gene_name", drop=TRUE]
sym <- as.character(sym)                  # avoid factor pitfalls
sym <- trimws(sym)
fallback_ids <- rownames(cts)
mask <- is.na(sym) | sym == "" | sym == "NA"
sym[mask] <- fallback_ids[mask]
sym <- trimws(sym)
mask2 <- is.na(sym) | sym == ""
sym[mask2] <- fallback_ids[mask2]
stopifnot(length(sym) == nrow(cts))
sym_unique <- make.unique(sym)            # character vector now
rownames(cts) <- sym_unique

options(Seurat.object.assay.version = "v5")
seurat_obj[["RNA_sym"]] <- CreateAssayObject(counts=cts)
DefaultAssay(seurat_obj) <- "RNA_sym"

# minimal preprocessing
seurat_obj <- NormalizeData(seurat_obj, normalization.method="LogNormalize", scale.factor=1e4, verbose=FALSE)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=3000, verbose=FALSE)
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj), verbose=FALSE)
seurat_obj <- RunPCA(seurat_obj, features=VariableFeatures(seurat_obj), npcs=30, verbose=FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims=1:30, k.param=30, verbose=FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution=0.4, algorithm=1, verbose=FALSE)
seurat_obj <- RunUMAP(seurat_obj, dims=1:30, n.neighbors=30, min.dist=0.3, verbose=FALSE)

p <- DimPlot(seurat_obj, group.by=groupby, label=TRUE) +
  ggtitle(paste0("UMAP — ", dataset, " / ", AGE_TAG, " / ", spefic_celltype)) + NoLegend()
ggsave(file.path(resultpath, paste0(PREFIX, "-UMAP_overview.png")), p, width=7, height=6, dpi=300, bg="white")

# ===== hdWGCNA =====
seurat_obj <- SetupForWGCNA(seurat_obj, gene_select="fraction", fraction=0.025, wgcna_name="tutorial")
seurat_obj <- MetacellsByGroups(seurat_obj, group.by=groupby, reduction="pca", k=kcells, max_shared=10, ident.group=groupby)
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(seurat_obj, group_name=spefic_celltype, group.by=groupby, assay="RNA_sym", layer="data")

seurat_obj <- TestSoftPowers(seurat_obj, networkType="signed")
plot_list <- PlotSoftPowers(seurat_obj); wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
write.csv(power_table, file.path(resultpath, paste0(PREFIX, "-softpower_table.csv")), row.names=FALSE)

seurat_obj <- ConstructNetwork(seurat_obj,
  tom_name=spefic_celltype, overwrite_tom=TRUE, soft_power=soft_power,
  deepSplit=deepSplit, detectCutHeight=detectCutHeight,
  minModuleSize=minModuleSize, mergeCutHeight=mergeCutHeight)

png(file.path(resultpath, paste0(PREFIX, "-dendrogram.png")), width=2400, height=1800, res=300)
PlotDendrogram(seurat_obj, main=paste0(PREFIX, " hdWGCNA Dendrogram"))
dev.off()

TOM <- GetTOM(seurat_obj)

seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj), verbose=FALSE)
seurat_obj <- ModuleEigengenes(seurat_obj)
hMEs <- GetMEs(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj, group.by=groupby, group_name=spefic_celltype)

seurat_obj <- ResetModuleNames(seurat_obj, new_name=paste0(dataset, "-", AGE_TAG, "-", spefic_celltype, "-M"))
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

  seurat_obj <- RunModuleUMAP(seurat_obj, n_hubs=10, n_neighbors=15, min_dist=0.1)
  png(file.path(resultpath, paste0(PREFIX, "-module-umap.png")), width=1800, height=1800, res=300)
  ModuleUMAPPlot(seurat_obj, edge.alpha=0.25, sample_edges=TRUE, edge_prop=0.1, label_hubs=2, keep_grey_edges=FALSE)
  dev.off()

  plot_module_native(obj=seurat_obj, resultpath=resultpath, dataset=dataset, age_tag=AGE_TAG,
                     spefic_celltype=spefic_celltype, group_var=groupby, n_genes=n_genes, assay="RNA_sym")

  module_enrich(modules_df=modules, result_path=file.path(resultpath, "enrich"),
                dataset=dataset, age_tag=AGE_TAG, spefic_celltype=spefic_celltype, cor_threshold=0.3)
} else {
  msg_sk("No non-grey modules detected; skip plotting/enrichment.")
}

msg_ok("All done → %s", resultpath)
