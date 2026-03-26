############################################################
##  CellAlign 对齐脚本  ——  读取共享 HVGs 列表（CSV 版）
############################################################
## ---------- ❶ 用户参数区 ---------------------------------
cds_a_path <- "D:/111/S-MG_subset_cds.rds"    # 样本 A
cds_b_path <- "D:/111/S-AG_subset_cds.rds"    # 样本 B
hvg_csv    <- "D:/111/HVGS-pairwise_intersections.csv"        # CSV 路径
hvg_col    <- "S-MG & S-AG"                     # 列名

celltypes_A <- c("MaSC","MaSC-t2-sg","LumSEC-Lac",
                 "LumSEC-Lip","Epi-Krt7")     # A 要保留的细胞类型
celltypes_B <- c("ASC-sg",'Epi-Lgals7','Epi-Krt7','LumSEC-AG-Pip')                          # B 全保留

celltype_col <- "newcelltype"                 # 细胞类型列名
ptime_col    <- "pseudotime"                  # 伪时序列名

win_sz   <- 0.1                               # 高斯窗口
num_pts  <- 200                               # 插值点数
############################################################

suppressPackageStartupMessages({
  library(monocle3)
  library(cellAlign)
  library(Matrix)
  library(data.table)      # 读 CSV 更快
})

## ---------- ❷ 读入 & 子集 -------------------------------
cdsA <- load_monocle_objects(cds_a_path)
cdsB <- load_monocle_objects(cds_b_path)

selA <- if ("all" %in% celltypes_A) rep(TRUE, ncol(cdsA)) else
  cdsA[[celltype_col]] %in% celltypes_A
selB <- if ("all" %in% celltypes_B) rep(TRUE, ncol(cdsB)) else
  cdsB[[celltype_col]] %in% celltypes_B

cdsA <- cdsA[, selA]
cdsB <- cdsB[, selB]

## ---------- ❸ 读取共享 HVGs ------------------------------
hvg_df   <- fread(hvg_csv, select = hvg_col)
genes_use <- unique(na.omit(hvg_df[[hvg_col]]))
genes_use <- intersect(genes_use,            # 再与真实存在基因取交集
                       intersect(rownames(cdsA), rownames(cdsB)))
message("Shared HVGs loaded: ", length(genes_use))

if (length(genes_use) < 50)
  stop("共享基因少于 50 个，检查 CSV 或列名是否正确。")

exprA <- as.matrix(exprs(cdsA)[genes_use, ])
exprB <- as.matrix(exprs(cdsB)[genes_use, ])

ptA  <- pseudotime(cdsA)
ptB  <- pseudotime(cdsB)
# === NEW: 清理伪时序向量 =================================

## ---------- ❹ 插值 & 缩放 -------------------------------
interA <- cellAlign::interWeights(exprA, ptA,
                                  winSz = win_sz, numPts = num_pts)
interB <- cellAlign::interWeights(exprB, ptB,
                                  winSz = win_sz, numPts = num_pts)

interA <- cellAlign::scaleInterpolate(interA)
interB <- cellAlign::scaleInterpolate(interB)

## ---------- ❺ 全局对齐 ---------------------------------
alignment <- cellAlign::globalAlign(interB$scaledData,interA$scaledData,
  scores  = list(query = interB$traj,
                 ref   = interA$traj),
  sigCalc = FALSE)

mapping <- cellAlign::mapRealDataGlobal(
  alignment,
  intTrajQuery = interB$traj, realTrajQuery = ptB,
  intTrajRef   = interA$traj, realTrajRef   = ptA)


## ---------- ❻ 写回 & 保存 -------------------------------
assign_aligned_pt <- function(cds, assign_list, meta_col, meta_table) {
  aligned_vec <- setNames(rep(NA_real_, ncol(cds)), colnames(cds))
  
  for (node_name in names(assign_list)) {
    cells_i <- assign_list[[node_name]]
    
    rows     <- which(meta_table[[meta_col]] == node_name)
    if (length(rows) == 0) next           # 理论上不会发生
    
    pts      <- unique(meta_table$ptRef[rows])
    if (length(pts) > 1) {                # 同名 meta‑node 不同 ptRef ⇒ 取平均
      warning(sprintf("meta‑node %s 有多个 ptRef (%.4f–%.4f)，取平均。",
                      node_name, min(pts), max(pts)))
      pts <- mean(pts)
    }
    aligned_vec[cells_i] <- pts[1]        # 同批细胞用同一坐标
  }
  colData(cds)$pseudotime_aligned <- aligned_vec
  cds
}

## --- 1.  Query (cdsB) ----------------------------------------
cdsB <- assign_aligned_pt(
  cds  = cdsB,
  assign_list = mapping$queryAssign,
  meta_col    = "metaNodeQuery",
  meta_table  = mapping$metaNodesPt)

## --- 2.  Reference (cdsA) ------------------------------------
cdsA <- assign_aligned_pt(
  cds  = cdsA,
  assign_list = mapping$refAssign,
  meta_col    = "metaNodeRef",
  meta_table  = mapping$metaNodesPt)


## -------- 保存 ------------------------------------------------
saveRDS(cdsA, sub("\\.rds$", "_cellalign.rds", cds_a_path))
saveRDS(cdsB, sub("\\.rds$", "_cellalign.rds", cds_b_path))
message("✓ pseudotime_aligned 已写入并保存 *_cellalign.rds")




###############################################################################
###############################################################################
##    Sliding‑window Spearman correlation: raw vs CellAlign‑aligned
###############################################################################
## ---------- 用户参数区 ------------------------------------------------------
raw_pt_col   <- "pseudotime"           # 未对齐伪时序列名
align_pt_col <- "pseudotime_aligned"   # CellAlign 对齐后伪时序列名

win_width <- 0.10      # 窗口宽度（归一化轴百分比）
step_size <- 0.01      # 步长
min_cells <- 5         # 每窗口最少细胞数（少于则设 NA）

out_prefix <- "MG_vs_AG_slidingCor"    # 输出图片前缀
###############################################################################

suppressPackageStartupMessages({
  library(Matrix)
  library(ggplot2)
  library(pbapply)      # 进度条
})

# ---------- 1. 准备工具函数 --------------------------------------------------
get_cor_curve <- function(exprA, exprB, ptA, ptB,
                          win = 0.1, step = 0.02, min_n = 5)
{
  ## 归一化到 0–1
  ptA <- (ptA - min(ptA, na.rm=TRUE)) / (max(ptA, na.rm=TRUE) - min(ptA, na.rm=TRUE))
  ptB <- (ptB - min(ptB, na.rm=TRUE)) / (max(ptB, na.rm=TRUE) - min(ptB, na.rm=TRUE))
  
  centers <- seq(0 + win/2, 1 - win/2, by = step)
  progress <- txtProgressBar(min = 0, max = length(centers), style = 3)
  rho <- numeric(length(centers))
  
  for (i in seq_along(centers)) {
    center <- centers[i]
    lo <- center - win/2
    hi <- center + win/2
    idxA <- which(ptA >= lo & ptA < hi)
    idxB <- which(ptB >= lo & ptB < hi)
    
    if (length(idxA) < min_n || length(idxB) < min_n) {
      rho[i] <- NA
    } else {
      vecA <- rowMeans(exprA[, idxA, drop = FALSE])
      vecB <- rowMeans(exprB[, idxB, drop = FALSE])
      rho[i] <- suppressWarnings(cor(vecA, vecB, method = "spearman"))
    }
    setTxtProgressBar(progress, i)
  }
  close(progress)
  data.frame(center = centers, rho = rho)
}

# ---------- 2. 提取表达矩阵 & 伪时序 ----------------------------------------
exprA <- as.matrix(exprs(cdsA))
exprB <- as.matrix(exprs(cdsB))

pt_raw_A   <- pseudotime(cdsA)
pt_raw_B   <- pseudotime(cdsB)
pt_align_A <- colData(cdsA)[[align_pt_col]]
pt_align_B <- colData(cdsB)[[align_pt_col]]

# ---------- 3. 计算滑窗相关 (raw & aligned) ---------------------------------
message("\n▶ 计算『未对齐伪时序』滑窗相关")
curve_raw <- get_cor_curve(exprA, exprB, pt_raw_A, pt_raw_B,
                           win_width, step_size, min_cells)

message("\n▶ 计算『CellAlign 对齐伪时序』滑窗相关")
curve_align <- get_cor_curve(exprA, exprB, pt_align_A, pt_align_B,
                             win_width, step_size, min_cells)

# ---------- 4. 画图 ---------------------------------------------------------
# 主题
p_theme <- theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))

## 图 1: raw
p1 <- ggplot(curve_raw, aes(center, rho)) +
  geom_line(color = "#d95f02", size = 1) +
  geom_point(color = "#d95f02", size = 1.5) +
  labs(title = "Sliding Spearman (raw pseudotime)",
       x = "Raw pseudotime (scaled 0–1)",
       y = "Spearman ρ") +
  ylim(0.8, 1) + p_theme

## 图 2: aligned
p2 <- ggplot(curve_align, aes(center, rho)) +
  geom_line(color = "#1b9e77", size = 1) +
  geom_point(color = "#1b9e77", size = 1.5) +
  labs(title = "Sliding Spearman (CellAlign‑aligned)",
       x = "Aligned pseudotime (0–1)",
       y = "Spearman ρ") +
  ylim(0.8, 1) + p_theme

## 图 3: overlay
overlay_df <- rbind(
  data.frame(curve_raw,   type = "raw"),
  data.frame(curve_align, type = "aligned")
)
p3 <- ggplot(overlay_df, aes(center, rho, color = type)) +
  geom_line(size = 1) + geom_point(size = 1.2) +
  scale_color_manual(values = c(raw = "#d95f02", aligned = "#1b9e77")) +
  labs(title = "Sliding Spearman: raw vs aligned",
       x = "Pseudotime (scaled)",
       y = "Spearman ρ",
       color = "") +
  ylim(0.85, 1) + p_theme

# ---------- 5. 输出 ---------------------------------------------------------
ggsave(paste0(out_prefix, "_raw.png"),  p1, width = 6, height = 4, dpi = 300,bg="white")
ggsave(paste0(out_prefix, "_aligned.png"),  p2, width = 6, height = 4, dpi = 300,bg="white")
ggsave(paste0(out_prefix, "_overlay.png"), p3, width = 6, height = 4, dpi = 300,bg="white")

message("\n✓ 三张图片已保存为 ",
        paste0(out_prefix, "_raw.png / _aligned.png / _overlay.png"))


