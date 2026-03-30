# ==============================================================================
# Script Name: plot_4quadrant_dynamic.R
# Description: 动态生成四象限富集图 (支持自定义路径、样本名自动组装)
# ==============================================================================

# ── 1. 加载依赖包 ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(scales) 
  library(grid)
})

# ── 2. 用户配置区 (请在此修改) ────────────────────────────────────────────────

# [配置 1] 设置通用根目录 (不要包含最后的 /)
root_dir <- "D:/Rfile/result/1225"

# [配置 2] 设置分析类型 (必须与文件名后缀一致)
# 例如文件名是 "ALL MG_GO_MF.csv"，这里就填 "GO_MF"
# 如果是 KEGG 分析，文件名通常是 "xxx_KEGG.csv"，这里填 "KEGG"
ANALYSIS_TYPE <- "KEGG" 

# [配置 3] 定义四个样本/分组的真实名称 (用于文件名拼接和图例显示)
# 对应关系：
# sample_1 -> 右上 (Top-Right)
# sample_2 -> 左上 (Top-Left)
# sample_3 -> 左下 (Bottom-Left)
# sample_4 -> 右下 (Bottom-Right)
sample_1 <- "ALL MG-purple"
sample_2 <- "ONLY 2AG"
sample_3 <- "ONLY 2AG&CG"
sample_4 <- "ONLY MG"

# [配置 4] 绘图参数
top_n_show  <- 5      # 每组标记多少个最显著的通路
point_alpha <- 0.8    # 点的透明度
base_size   <- 14     # 基础字体大小

# [配置 5] 定义颜色 (顺序对应 sample_1 到 sample_4)
# 依次为: 右上(红), 左上(蓝), 左下(绿), 右下(紫)
colors_vec <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488")

# ──────────────────────────────────────────────────────────────────────────────
# 自动逻辑区 (以下部分通常不需要修改)
# ──────────────────────────────────────────────────────────────────────────────

# 1. 自动组装输入文件路径
# 假设文件结构为: root_dir/tables/[SampleName]_[AnalysisType].csv
files_list <- list()
files_list[[sample_1]] <- file.path(root_dir, "tables", paste0(sample_1, "_", ANALYSIS_TYPE, ".csv"))
files_list[[sample_2]] <- file.path(root_dir, "tables", paste0(sample_2, "_", ANALYSIS_TYPE, ".csv"))
files_list[[sample_3]] <- file.path(root_dir, "tables", paste0(sample_3, "_", ANALYSIS_TYPE, ".csv"))
files_list[[sample_4]] <- file.path(root_dir, "tables", paste0(sample_4, "_", ANALYSIS_TYPE, ".csv"))

# 2. 自动设置输出目录和文件名
out_dir <- file.path(root_dir, "plots_4quadrant")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# 输出文件名：包含分析类型
out_name <- paste0("4Quad_Enrichment_purple", ANALYSIS_TYPE, ".pdf")
out_path <- file.path(out_dir, out_name)

# 3. 颜色名映射
names(colors_vec) <- c(sample_1, sample_2, sample_3, sample_4)

# ── 3. 数据处理函数 ───────────────────────────────────────────────────────────

parse_ratio <- function(x) {
  if(is.numeric(x)) return(x)
  if(all(is.na(x))) return(0)
  sapply(x, function(val) {
    if(grepl("/", val)) {
      parts <- as.numeric(str_split(val, "/", simplify = TRUE))
      return(parts[1] / parts[2])
    } else {
      return(as.numeric(val))
    }
  })
}

load_enrich_data <- function(filepath, group_name) {
  if (!file.exists(filepath)) {
    warning(paste("⚠️ 文件不存在，跳过:", filepath))
    return(NULL)
  }
  
  df <- read_csv(filepath, show_col_types = FALSE)
  
  # 兼容 p.adjust / qvalue / PValue
  if ("p.adjust" %in% names(df)) {
    df$P_Use <- df$p.adjust
  } else if ("qvalue" %in% names(df)) {
    df$P_Use <- df$qvalue
  } else if ("PValue" %in% names(df)) {
    df$P_Use <- df$PValue
  } else {
    stop(paste("文件", filepath, "中未找到 p.adjust/qvalue/PValue 列"))
  }
  
  df %>%
    mutate(
      Group = group_name,
      Description = str_remove(Description, " - Mus musculus \\(house mouse\\)"),
      GeneRatio_Num = parse_ratio(GeneRatio),
      logP = -log10(P_Use)
    ) %>%
    filter(P_Use < 0.05) %>% 
    select(ID, Description, Group, GeneRatio_Num, logP, Count, P_Use)
}

# ── 4. 数据加载与象限映射 ─────────────────────────────────────────────────────

message(">>> 正在加载数据...")
all_data <- list()
quadrant_order <- names(files_list) # 这里现在是真实的样本名

for (grp in quadrant_order) {
  message(paste("   Reading:", grp))
  all_data[[grp]] <- load_enrich_data(files_list[[grp]], grp)
}
combined_df <- bind_rows(all_data)

# 锁定因子顺序，保证图例顺序一致
combined_df$Group <- factor(combined_df$Group, levels = quadrant_order)

message(">>> 计算象限坐标...")
# 逻辑: 1=右上(++), 2=左上(-+), 3=左下(--), 4=右下(+-)
plot_df <- combined_df %>%
  mutate(
    x_plot = case_when(
      Group == quadrant_order[1] ~ logP,
      Group == quadrant_order[2] ~ -logP,
      Group == quadrant_order[3] ~ -logP,
      Group == quadrant_order[4] ~ logP
    ),
    y_plot = case_when(
      Group == quadrant_order[1] ~ GeneRatio_Num,
      Group == quadrant_order[2] ~ GeneRatio_Num,
      Group == quadrant_order[3] ~ -GeneRatio_Num,
      Group == quadrant_order[4] ~ -GeneRatio_Num
    )
  )

# 提取 Top N
label_df <- plot_df %>%
  group_by(Group) %>%
  slice_min(P_Use, n = top_n_show) %>%
  ungroup()

# 计算坐标轴最大值 (对称)
max_x <- max(abs(plot_df$x_plot), na.rm=TRUE) * 1.1
max_y <- max(abs(plot_df$y_plot), na.rm=TRUE) * 1.1

# ── 5. 绘图 ───────────────────────────────────────────────────────────────────

message(">>> 开始绘图...")

# 动态标题
plot_title <- paste0("Enrichment Map: ", ANALYSIS_TYPE)

p <- ggplot(plot_df, aes(x = x_plot, y = y_plot)) +
  
  # 1. 坐标参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  
  # 2. 散点 (使用真实样本名)
  geom_point(aes(color = Group, size = Count), alpha = point_alpha) +
  
  # 3. 标签
  geom_text_repel(data = label_df, 
                  aes(label = Description), 
                  color = "black",          
                  size = 3.5, 
                  fontface = "italic",
                  bg.color = "white", bg.r = 0.15,
                  box.padding = 0.5,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  
  # 4. 坐标轴 (隐藏负号)
  scale_x_continuous(labels = abs, limits = c(-max_x, max_x),
                     name = paste0("-log10(p.adjust)  [", ANALYSIS_TYPE, "]")) +
  scale_y_continuous(labels = abs, limits = c(-max_y, max_y),
                     name = "Gene Ratio") +
  
  # 5. 颜色与图例 (使用真实样本名)
  scale_color_manual(
    values = colors_vec,
    name = "Samples" # 图例标题
  ) +
  scale_size_continuous(
    range = c(2, 7), 
    name = "Gene Count"
  ) +
  
  # 6. 主题
  theme_bw(base_size = base_size) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  # 图例点加大
  guides(color = guide_legend(override.aes = list(size = 5))) +
  
  labs(title = plot_title)

# ── 6. 保存 ───────────────────────────────────────────────────────────────────

ggsave(out_path, p, width = 11, height = 9)

message(paste("✅ 完成！"))
message(paste("   分析类型:", ANALYSIS_TYPE))
message(paste("   图片路径:", out_path))