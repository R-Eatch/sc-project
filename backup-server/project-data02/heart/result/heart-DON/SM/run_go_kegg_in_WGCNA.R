#!/usr/bin/env Rscript

# =======================
# 全局配置（请按需修改）
# =======================
dataset        <- "DON"                                    # 数据集标签
celltype       <- "Smooth_Muscle"                          # 细胞类型标签
modules_csv    <- paste0("./result",dataset,"-",celltype,"modules.csv")  # GetModules 导出的模块表
outdir         <- "./result/enrich"                        # 输出根目录（将含 plots/ 与 tables/）
organism       <- "human"                                  # "human" 或 "mouse"
cor_threshold  <- 0.3                                      # |kME| 选择阈值
min_gene_count <- 3                                        # 最少基因数门槛
top_show       <- 15                                       # dotplot显示的最多条目
include_kegg   <- FALSE                                     # ← 新增：是否执行 KEGG
plot_width     <- 14
plot_height    <- 14
gene_col_opt   <- NULL                                     # 可设为 "gene_name" 或 "gene"；NULL 自动检测
kme_prefix     <- "kME_"                                   # 与 GetModules 输出一致

# =======================
# 依赖加载
# =======================
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(clusterProfiler)
  library(enrichplot)   # dotplot 可视化
})

# =======================
# 工具函数
# =======================
empty_plot <- function(txt) {
  ggplot() +
    annotate("text", x = .5, y = .5, label = txt, size = 4, lineheight = .9) +
    theme_void()
}

save_tbl <- function(obj, path) write.csv(as.data.frame(obj), path, row.names = FALSE)

# 如果 cairo_pdf 不可用，自动回退
ggsave_pdf_safe <- function(filename, plot, width, height) {
  ok <- TRUE
  tryCatch(
    ggsave(filename, plot, width = width, height = height, device = cairo_pdf),
    error = function(e) {
      message("⚠️  cairo_pdf 不可用，使用默认 pdf 设备：", e$message)
      ok <<- FALSE
    }
  )
  if (!ok) {
    alt <- gsub("\\.pdf$", "_fallback.pdf", filename)
    ggsave(alt, plot, width = width, height = height)
  }
}

# =======================
# 主富集函数（对 data.frame 执行）
# =======================
module_enrich_from_df <- function(
    modules_df,
    result_path,
    dataset,
    spefic_celltype,
    organism          = c("mouse","human"),
    gene_col          = NULL,
    kme_prefix        = "kME_",
    cor_threshold     = 0.5,
    min_gene_count    = 3,
    top_show          = 15,
    include_kegg      = TRUE,
    go_p_cutoff       = 0.05,
    go_q_cutoff       = 0.20,
    kegg_p_cutoff     = 0.10,
    kegg_q_cutoff     = 0.25,
    plot_width        = 14,
    plot_height       = 14
) {
  suppressPackageStartupMessages({
    library(dplyr); library(ggplot2); library(patchwork)
    library(clusterProfiler)
  })
  
  # ---------- 物种 ----------
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
  
  # ---------- 列名 ----------
  if (is.null(gene_col)) {
    if ("gene_name" %in% names(modules_df)) gene_col <- "gene_name"
    else if ("gene" %in% names(modules_df)) gene_col <- "gene"
    else stop("找不到基因名列，请设置 gene_col 或确保存在 gene_name / gene 列。")
  }
  if (!"module" %in% names(modules_df)) stop("modules_df 缺少列：module")
  
  # ---------- 输出目录 ----------
  plot_dir  <- file.path(result_path, "plots")
  table_dir <- file.path(result_path, "tables")
  dir.create(plot_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ---------- 占位与兜底 ----------
  empty_plot <- function(txt) {
    ggplot() +
      annotate("text", x = .5, y = .5, label = txt, size = 4, lineheight = .9) +
      coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +  # 明确坐标，避免 scale 计算拿到 NULL
      theme_void()
  }
  force_plot <- function(p, msg = "empty") {
    if (inherits(p, "ggplot") || inherits(p, "patchwork")) return(p)
    empty_plot(msg)
  }
  save_tbl <- function(obj, path) write.csv(as.data.frame(obj), path, row.names = FALSE)
  
  # ---------- 遍历模块 ----------
  all_modules <- unique(modules_df$module)
  message("🟢 共检测到 ", length(all_modules), " 个模块：", paste(all_modules, collapse = ", "))
  
  for (mod in all_modules) {
    kme_col <- paste0(kme_prefix, mod)
    this_df <- modules_df %>% dplyr::filter(.data$module == !!mod)
    
    if (!kme_col %in% names(this_df)) {
      warning("模块 ", mod, ": 找不到 kME 列 ", kme_col, "，跳过。")
      next
    }
    
    genes <- this_df %>%
      dplyr::filter(abs(.data[[kme_col]]) >= cor_threshold) %>%
      dplyr::pull(.data[[gene_col]]) %>%
      unique() %>% na.omit()
    
    if (length(genes) < min_gene_count) {
      warning(sprintf("模块 %s: 基因数(%d) < min_gene_count(%d)，跳过。", mod, length(genes), min_gene_count))
      next
    }
    message(sprintf("  • 模块 %s：用于富集的基因数 = %d", mod, length(genes)))
    
    # ---- 1) GO（先放占位，后覆盖）----
    go_plots <- setNames(
      list(empty_plot(paste(mod, "no GO BP")),
           empty_plot(paste(mod, "no GO CC")),
           empty_plot(paste(mod, "no GO MF"))),
      c("BP","CC","MF")
    )
    
    for (ont in c("BP","CC","MF")) {
      go_res <- tryCatch(
        enrichGO(
          gene          = genes,
          OrgDb         = OrgDb,
          keyType       = "SYMBOL",
          ont           = ont,
          pAdjustMethod = "BH",
          pvalueCutoff  = go_p_cutoff,
          qvalueCutoff  = go_q_cutoff
        ),
        error = function(e) {
          warning("模块 ", mod, " GO-", ont, " 失败: ", e$message); return(NULL)
        }
      )
      
      if (!is.null(go_res) && nrow(go_res@result) > 0) {
        save_tbl(go_res@result, file.path(table_dir, sprintf("%s-%s-%s_GO_%s.csv", dataset, spefic_celltype, mod, ont)))
        go_plots[[ont]] <- tryCatch(
          enrichplot::dotplot(go_res, showCategory = top_show) +
            ggtitle(sprintf("%s-%s %s – GO-%s", dataset, spefic_celltype, mod, ont)) +
            theme_classic(base_size = 11),
          error = function(e) {
            warning("模块 ", mod, " GO-", ont, " 绘图失败: ", e$message)
            empty_plot(paste(mod, "GO-", ont, " plot error"))
          }
        )
      }
    }
    
    # ---- 2) KEGG（可选，先放占位）----
    kegg_plot <- empty_plot(paste(mod, "KEGG disabled"))
    if (isTRUE(include_kegg)) {
      entrez_map <- suppressMessages(bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb))
      if (!is.null(entrez_map) && nrow(entrez_map) > 0) {
        entrez_ids <- unique(entrez_map$ENTREZID)
        kegg_res <- tryCatch(
          enrichKEGG(
            gene          = entrez_ids,
            organism      = kegg_org,
            keyType       = "ncbi-geneid",
            pAdjustMethod = "BH",
            pvalueCutoff  = kegg_p_cutoff,
            qvalueCutoff  = kegg_q_cutoff
          ),
          error = function(e) { warning("模块 ", mod, " KEGG 失败: ", e$message); NULL }
        )
        
        if (!is.null(kegg_res) && methods::is(kegg_res, "enrichResult")) {
          kegg_df <- tryCatch(as.data.frame(kegg_res), error = function(e) NULL)
          if (!is.null(kegg_df) && nrow(kegg_df) > 0) {
            kegg_df$geneID <- vapply(kegg_df$geneID, function(x) {
              ids  <- strsplit(x, "/")[[1]]
              syms <- entrez_map$SYMBOL[match(ids, entrez_map$ENTREZID)]
              paste(na.omit(syms), collapse = "/")
            }, character(1))
            save_tbl(kegg_df, file.path(table_dir, sprintf("%s-%s-%s_KEGG.csv", dataset, spefic_celltype, mod)))
            
            kegg_plot <- tryCatch(
              enrichplot::dotplot(kegg_res, showCategory = top_show) +
                ggtitle(sprintf("%s-%s %s – KEGG", dataset, spefic_celltype, mod)) +
                theme_classic(base_size = 11),
              error = function(e) {
                warning("模块 ", mod, " KEGG 绘图失败: ", e$message)
                empty_plot(paste(mod, "KEGG plot error"))
              }
            )
          } else {
            kegg_plot <- empty_plot(paste(mod, "no significant KEGG"))
          }
        }
      } else {
        kegg_plot <- empty_plot(paste(mod, "no valid genes for KEGG"))
      }
    }
    
    # ---- 3) 组合与保存（统一兜底）----
    combined <- patchwork::wrap_plots(
      list(
        force_plot(go_plots$BP, "GO-BP empty"),
        force_plot(go_plots$CC, "GO-CC empty"),
        force_plot(go_plots$MF, "GO-MF empty"),
        force_plot(kegg_plot,   "KEGG empty")
      ),
      ncol = 2
    )
    
    pdf_path <- file.path(plot_dir, sprintf("%s-%s-%s_GO_KEGG.pdf", dataset, spefic_celltype, mod))
    tryCatch(
      ggsave(pdf_path, combined, width = plot_width, height = plot_height, device = cairo_pdf),
      error = function(e) {
        warning("保存 PDF 失败（cairo）：", e$message, "，尝试默认 pdf 设备。")
        tryCatch(
          ggsave(gsub("\\.pdf$", "_fallback.pdf", pdf_path), combined, width = plot_width, height = plot_height),
          error = function(e2) warning("保存 PDF 仍失败：", e2$message)
        )
      }
    )
  }
  
  message("全部完成！结果已保存至：", result_path)
  invisible(TRUE)
}


# =======================
# 读取 CSV 并运行
# =======================
if (!file.exists(modules_csv)) stop("找不到文件：", modules_csv)
modules_df <- suppressMessages(readr::read_csv(modules_csv, show_col_types = FALSE))

module_enrich_from_df(
  modules_df       = modules_df,
  result_path      = outdir,
  dataset          = dataset,
  spefic_celltype  = celltype,
  organism         = organism,
  gene_col         = gene_col_opt,
  kme_prefix       = kme_prefix,
  cor_threshold    = cor_threshold,
  min_gene_count   = min_gene_count,
  top_show         = top_show,
  include_kegg     = include_kegg,
  plot_width       = plot_width,
  plot_height      = plot_height
)
