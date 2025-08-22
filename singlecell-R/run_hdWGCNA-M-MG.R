# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
##############################################
## global variable
groupby <- 'newcelltype'
harmony_vars <- 'sample'
reduction <- 'pca_harmony'
spefic_celltype <- 'Lum-Kit'
resultpath <- 'D:/Data/hdwgcna-result/human-heart/FB/'
RDSpath <- 'D:/Data/NormalHeart_cleaned.rds'
dataset <- 'M-MG'
n_genes <- 5
nThreads <- 4
soft_power <- NULL
deepSplit       = 2    
detectCutHeight = 0.999
minModuleSize   = 100
mergeCutHeight  = 0.15
##############################################
## settings
dir.create(resultpath, recursive = TRUE, showWarnings = FALSE)
# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = nThreads)
##############################################
## function

plot_module_native <- function(obj,
                               resultpath,
                               dataset,
                               spefic_celltype,
                               n_genes      = 5,
                               group_var  = "newcelltype",
                               reduction  = "umap",
                               assay      = "RNA") {
  ## ── 依赖包 ────────────────────────────────────────────────
  suppressPackageStartupMessages({
    library(dplyr)
    library(Seurat)
    library(patchwork)
    library(ggplot2)
  })
  
  ## ── 输出目录 ─────────────────────────────────────────────
  if (!dir.exists(resultpath)) dir.create(resultpath, recursive = TRUE)
  if (!grepl("/$", resultpath)) resultpath <- paste0(resultpath, "/")
  obj <- seurat_obj
  ## ── 1. 拿模块表 & 取 hub genes ───────────────────────────
  modules_tbl <- GetModules(obj) %>%
    dplyr::filter(module != "grey") %>%
    dplyr::select(gene_name, module, dplyr::starts_with("kME_"))
  
  hub_tbl <- modules_tbl %>% 
    dplyr::group_by(module) %>% 
    dplyr::group_modify(~{
      col_kme <- paste0("kME_", .y$module)
      .x %>% 
        dplyr::arrange(dplyr::desc(.data[[col_kme]])) %>% 
        dplyr::slice_head(n = n_genes)   # ← 关键：.env 取函数环境变量
    }) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(gene_tag = paste0(gene_name, " (", module, ")"))
  fp_list <- list()
  mod_list <- split(hub_tbl, hub_tbl$module)
  for (mod in names(mod_list)) {
    df_mod    <- mod_list[[mod]]
    genes_mod <- df_mod$gene_name
    
    
    if (length(genes_mod) == 0) {
      warning("模块 ", mod, " 无可用基因，跳过")
      next
    }
    
    message("绘制模块 ", mod, ": ",
            paste(genes_mod, collapse = ", "))
    
    fp_list[[mod]] <- FeaturePlot(
      obj,
      features  = genes_mod,
      reduction = reduction,
      ncol      = length(genes_mod),
      order     = TRUE
    ) +
      plot_annotation(title = paste0(mod, " hub genes")) +
      theme(plot.title = element_text(hjust = .5))
  }
  p_feature <- patchwork::wrap_plots(fp_list, ncol = 1)
  
  ## ── 3. DotPlot（加模块标签）──────────────────────────────
  gene_order    <- hub_tbl$gene_name                 # 保持顺序
  gene_labels   <- hub_tbl$gene_tag
  names(gene_labels) <- gene_order                   # 供 scale_y
  
  p_dot <- DotPlot(obj,
                   features = gene_order,
                   group.by = group_var,
                   assay    = assay) +
    scale_y_discrete(labels = gene_labels) +
    labs(y = "Gene (module)", x = group_var,
         title = "DotPlot of hub genes") +
    theme(axis.text.y = element_text(size = 8))
  
  ## ── 4. 保存 ──────────────────────────────────────────────
  file_fp  <- paste0(resultpath, dataset, "-", spefic_celltype, "-FeaturePlot.png")
  file_dot <- paste0(resultpath, dataset, "-", spefic_celltype, "-DotPlot.png")
  
  ggsave(file_fp,  p_feature, width = 15, height = 2 + 2.5 * length(fp_list), dpi = 300,bg = 'white',limitsize = FALSE)
  ggsave(file_dot, p_dot,     width = 10 + 2.5 * length(fp_list), height = 6,                dpi = 300,bg = 'white',limitsize = FALSE)
  
  message("✅ Saved:\n  ", file_fp, "\n  ", file_dot)
  invisible(list(feature = p_feature,
                 dot     = p_dot,
                 hub_genes = hub_tbl))
}

# ==============================================================
# Function: module_enrich
# 对 hdWGCNA "modules" 表的每个模块做 GO/KEGG 富集
# ==============================================================

module_enrich <- function(
    modules_df,                # ← GetModules(seurat_obj) %>% subset(module!="grey")
    result_path,               # ← 保存根目录，如 "results/enrich/"
    dataset,                   # ← 数据集标签，如 "S-MG"
    spefic_celltype,           # ← 细胞类型标签，如 "Epi-Krt7"
    organism       = c("mouse", "human"), # 用哪个注释库
    gene_col       = "gene_name",         # 基因名列
    cor_threshold  = 0.5,                 # |kME| 阈值
    top_show       = 15,                  # dotplot 最多展示多少条目
    plot_width     = 14,                  # PDF 尺寸
    plot_height    = 14
) {
  # ---------- 0. 依赖包 ----------
  suppressPackageStartupMessages({
    library(dplyr)
    library(clusterProfiler)
    library(ggplot2)
    library(patchwork)
  })
  
  # ---------- 1. 物种信息 ----------
  organism <- match.arg(organism)
  if (organism == "mouse") {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    OrgDb    <- org.Mm.eg.db
    kegg_org <- "mmu"
  } else {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    OrgDb    <- org.Hs.eg.db
    kegg_org <- "hsa"
  }
  
  # ---------- 2. 输出目录 ----------
  plot_dir  <- file.path(result_path, "plots")
  table_dir <- file.path(result_path, "tables")
  dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- 3. 辅助：空图 / 写表 ----------
  empty_plot <- function(txt) {
    ggplot() + annotate("text", x = .5, y = .5, label = txt, size = 4) +
      theme_void()
  }
  save_tbl <- function(obj, path) write.csv(as.data.frame(obj), path, row.names = FALSE)
  
  # ---------- 4. 遍历模块 ----------
  all_modules <- unique(modules_df$module)
  message("🟢  共检测到 ", length(all_modules), " 个模块")
  
  for (mod in all_modules) {
    kme_col <- paste0("kME_", mod)
    this_df <- modules_df %>% filter(module == !!mod)
    
    if (!kme_col %in% names(this_df)) {
      warning("找不到列 ", kme_col, "，跳过模块 ", mod)
      next
    }
    
    # 4.1 过滤基因
    genes <- this_df %>% 
      filter(abs(.data[[kme_col]]) >= cor_threshold) %>% 
      pull(.data[[gene_col]]) %>% 
      unique() %>% 
      na.omit()
    
    if (length(genes) == 0) {
      warning("模块 ", mod, " 在阈值 ", cor_threshold,
              " 下无基因满足条件，跳过。")
      next
    }
    
    message("  • 模块 ", mod, " → ", length(genes), " genes")
    
    # 4.2 占位 plot 容器
    go_plots <- vector("list", 3); names(go_plots) <- c("BP","CC","MF")
    
    # 4.3 GO 富集
    for (ont in c("BP","CC","MF")) {
      go_res <- tryCatch(
        enrichGO(gene = genes, OrgDb = OrgDb, keyType = "SYMBOL",
                 ont = ont, pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.20),
        error = function(e) NULL
      )
      
      if (!is.null(go_res) && nrow(go_res)) {
        save_tbl(go_res@result,
                 file.path(table_dir,
                           sprintf("%s-%s-%s_GO_%s.csv",
                                   dataset, spefic_celltype, mod, ont)))
        
        go_res@result <- head(go_res[order(go_res$p.adjust), ], top_show)
        go_plots[[ont]] <- dotplot(go_res, showCategory = top_show) +
          ggtitle(sprintf("%s-%s %s – GO-%s",
                          dataset, spefic_celltype, mod, ont)) +
          theme_classic(base_size = 11)
      } else {
        go_plots[[ont]] <- empty_plot(paste(mod, "no GO", ont))
      }
    }
    
    # 4.4 KEGG 富集
    entrez_map <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                       OrgDb = OrgDb)
    entrez_ids <- unique(entrez_map$ENTREZID)
    
    kegg_res <- tryCatch(
      enrichKEGG(gene = entrez_ids, organism = kegg_org,
                 keyType = "ncbi-geneid",
                 pAdjustMethod = "BH", pvalueCutoff = 0.10),
      error = function(e) NULL
    )
    
    if (!is.null(kegg_res) && nrow(kegg_res)) {
      # 把 EntrezID 转回 Symbol
      kegg_res@result$geneID <- sapply(kegg_res@result$geneID, function(x) {
        ids  <- strsplit(x, "/")[[1]]
        syms <- entrez_map$SYMBOL[match(ids, entrez_map$ENTREZID)]
        paste(syms, collapse = "/")
      })
      save_tbl(kegg_res@result,
               file.path(table_dir,
                         sprintf("%s-%s-%s_KEGG.csv",
                                 dataset, spefic_celltype, mod)))
      
      kegg_res@result <- head(kegg_res[order(kegg_res$p.adjust), ], top_show)
      kegg_plot <- dotplot(kegg_res, showCategory = top_show) +
        ggtitle(sprintf("%s-%s %s – KEGG",
                        dataset, spefic_celltype, mod)) +
        theme_classic(base_size = 11)
    } else {
      kegg_plot <- empty_plot(paste(mod, "no KEGG"))
    }
    
    # 4.5 组合并保存 PDF
    combined <- (go_plots[["BP"]] | go_plots[["CC"]]) /
      (go_plots[["MF"]] | kegg_plot)
    
    pdf_path <- file.path(
      plot_dir,
      sprintf("%s-%s-%s_GO_KEGG.pdf",
              dataset, spefic_celltype, mod)
    )
    ggsave(pdf_path, combined,
           width = plot_width, height = plot_height, device = cairo_pdf)
  }
  
  cat("🎉  module_enrich 全部完成！结果已保存至：", result_path, "\n")
  invisible(TRUE)
}

##############################################
## hdWGCNA
seurat_obj <- readRDS(RDSpath)
seurat_obj <- SeuratObject::UpdateSeuratObject(seurat_obj)

p <- DimPlot(seurat_obj, group.by=groupby, label=TRUE) +
  umap_theme() + ggtitle('celltype') + NoLegend()
ggsave(filename =paste0(resultpath,dataset,'-',spefic_celltype,'-celltype.png') ,plot = p)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = groupby, # specify the columns in seurat_obj@meta.data to group by
  reduction = reduction, # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = groupby # set the Idents of the metacell seurat object
)

seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = spefic_celltype, # the name of the group of interest in the group.by column
  group.by=groupby, # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)
# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)
# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)
# assemble with patchwork
wrap_plots(plot_list, ncol=2)
power_table <- GetPowerTable(seurat_obj)
head(power_table)

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = spefic_celltype,overwrite_tom = TRUE,soft_power =  soft_power,
  deepSplit       = deepSplit,
  detectCutHeight = detectCutHeight,  
  minModuleSize   = minModuleSize,     
  mergeCutHeight  = mergeCutHeight,   
)
png(paste0(resultpath,dataset,'-',spefic_celltype,"dendrogram.png"), width = 2400, height = 1800, res = 300)
PlotDendrogram(seurat_obj, main=paste0(dataset, '-',spefic_celltype,' ','hdWGCNA Dendrogram'))
dev.off()
TOM <- GetTOM(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars=harmony_vars
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)
# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by =groupby, group_name = spefic_celltype
)
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(dataset,'-',spefic_celltype,'-M')
)
saveRDS(seurat_obj, file=paste0(resultpath,dataset,'-',spefic_celltype,'-hdWGCNA_object.rds'))

#######################visualization###################
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')
p2 <- PlotKMEs(seurat_obj, ncol=5,n_hubs = 15)
ggsave(filename =paste0(resultpath,dataset,'-',spefic_celltype,'KMEsplot.png') ,plot = p2,dpi = 300,bg = 'white',width = 30,height = 15)
write.csv(modules,file = paste0(resultpath,dataset,'-',spefic_celltype,'modules.csv'),row.names = FALSE)

ModuleNetworkPlot(
  seurat_obj,
  outdir = paste0(resultpath,'moduleNetwork-default')
)
ModuleNetworkPlot(
  seurat_obj, 
  outdir=paste0(resultpath,'moduleNetwork-20'), # new folder name
  n_inner = 20, # number of genes in inner ring
  n_outer = 30, # number of genes in outer ring
  n_conns = Inf, # show all of the connections
  plot_size=c(10,10), # larger plotting area
  vertex.label.cex=1 # font size
)

seurat_obj <- readRDS(paste0(resultpath,dataset,'-',spefic_celltype,'-hdWGCNA_object.rds'))
# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)
p3 <- wrap_plots(plot_list, ncol=5)
ggsave(filename =paste0(resultpath,dataset,'-',spefic_celltype,'hMEs-fp.png') ,plot = p3,dpi = 200,bg = 'white',width = 30,height = length(modules)*1,limitsize = FALSE)

modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']
png(paste0(resultpath,dataset,'-',spefic_celltype,"HubGeneNetwork.png"), width = 2400, height = 1800, res = 300)
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods
)
dev.off()
# modules umap
if(FALSE){
  seurat_obj <- RunModuleUMAP(
    seurat_obj,
    n_hubs = 10, # number of hub genes to include for the UMAP embedding
    n_neighbors=15, # neighbors parameter for UMAP
    min_dist=0.1 # min distance between points in UMAP space
  )
  png(paste0(resultpath,dataset,'-',spefic_celltype,"module-umap.png"), width = 1800, height = 1800, res = 300)
  ModuleUMAPPlot(
    seurat_obj,
    edge.alpha=0.25,
    sample_edges=TRUE,
    edge_prop=0.1, # proportion of edges to sample (20% here)
    label_hubs=2 ,# how many hub genes to plot per module?
    keep_grey_edges=FALSE
  )
  dev.off()
}
#####################################################
seurat_obj <- readRDS(paste0(resultpath,dataset,'-',spefic_celltype,'-hdWGCNA_object.rds'))
#####################################################
## plot featureplot
plot_module_native(
  obj = seurat_obj,
  resultpath      = resultpath,   # 会自动创建
  dataset         = dataset,
  spefic_celltype = spefic_celltype,
  group_var       = groupby,
  n_genes= n_genes
)
#####################################################
modules <- GetModules(seurat_obj) %>%          # hdWGCNA 输出
  subset(module != "grey")            # 去掉 grey 模块
## GO KEGG enrich
module_enrich(
  modules_df      = modules,
  result_path     = paste0(resultpath,"enrich/"),
  dataset         = dataset,
  spefic_celltype = spefic_celltype,
  organism        = "mouse",     # 或 "human"
  cor_threshold   = 0.5          # 只保留 |kME| ≥ 0.6 的 hub genes
)
#####################################################


