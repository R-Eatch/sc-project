# --- 1. 加载必要的库 ---
# 确保已安装这些包: install.packages(c("readr", "dplyr", "tidyr", "ggplot2", "ggrepel", "RColorBrewer"))
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
dir_name <- 'D:/111/DEG_Analysis_Output/GO-KEGG-OUTPUT/'
dataset <- 'juvenile'
# --- 2. 读取数据 ---
go_file <- paste0(dir_name,"GO-result-ALL_MG-juvenile.csv")
kegg_file <- paste0(dir_name,"KEGG-result-ALL_MG-juvenile.csv")

if (!file.exists(go_file)) stop("错误：找不到 GO 结果文件: ", go_file)
if (!file.exists(kegg_file)) stop("错误：找不到 KEGG 结果文件: ", kegg_file)

go_data <- read_csv(go_file)
kegg_data <- read_csv(kegg_file)
kegg_data$Description <- gsub(
  pattern     = " - Mus musculus \\(house mouse\\)",
  replacement = "",
  x           = kegg_data$Description
)
# --- 3. 定义需要显示/高亮的 GO 和 KEGG ID ---
# !!! 关键步骤：请根据您的分析结果修改这些 ID 列表 !!!
# 这个列表同时用于高亮显示点（两个版本）和添加标签（仅标签版本）
go_ids_to_select <- c(
  "GO:0022613",  # ribonucleoprotein complex biogenesis
  "GO:0042254",  # ribosome biogenesis
  "GO:0006260",  # DNA replication
  "GO:0034470",  # ncRNA processing
  "GO:0008380"   # RNA splicing
)

kegg_ids_to_select <- c(
  "mmu03030",  # DNA replication
  "mmu03040",  # Spliceosome
  "mmu04110",  # Cell cycle
  "mmu03050",  # Proteasome
  "mmu03013"   # RNA transport
)

# --- 4. 数据预处理 ---
# (函数与之前相同)
process_enrichment_data <- function(df, type, x_axis_sign = 1) {
  if (!"GeneRatio" %in% names(df)) stop("错误：数据框中缺少 'GeneRatio' 列。")
  if (!"qvalue" %in% names(df)) stop("错误：数据框中缺少 'qvalue' 列。")
  if (!"ID" %in% names(df)) stop("错误：数据框中缺少 'ID' 列。")
  if (!"Description" %in% names(df)) stop("错误：数据框中缺少 'Description' 列。")
  if (!"Count" %in% names(df)) stop("错误：数据框中缺少 'Count' 列。")
  
  df %>%
    select(ID, Description, GeneRatio, qvalue, Count) %>%
    filter(grepl("^\\d+/\\d+$", GeneRatio)) %>%
    separate(GeneRatio, into = c("numerator", "denominator"), sep = "/", convert = TRUE, remove = FALSE) %>%
    mutate(GeneRatioNumeric = numerator / denominator) %>%
    mutate(qvalue = if_else(qvalue == 0 | is.na(qvalue), 1e-300, qvalue)) %>%
    mutate(negLog10q = -log10(qvalue)) %>%
    mutate(x_axis_val = negLog10q * x_axis_sign) %>%
    mutate(Type = type) %>%
    mutate(Description = as.character(Description)) %>%
    select(ID, Description, GeneRatioNumeric, qvalue, negLog10q, x_axis_val, Count, Type)
}

go_processed <- process_enrichment_data(go_data, type = "GO", x_axis_sign = 1)
kegg_processed <- process_enrichment_data(kegg_data, type = "KEGG", x_axis_sign = -1)

combined_data <- bind_rows(go_processed, kegg_processed) %>%
  filter(!is.na(x_axis_val), is.finite(x_axis_val), !is.na(GeneRatioNumeric))

# --- 5. 准备高亮/标签数据 ---

# 筛选需要高亮显示的点的数据 (两个版本都需要)
highlighting_data <- combined_data %>%
  filter((Type == "GO" & ID %in% go_ids_to_select) |
           (Type == "KEGG" & ID %in% kegg_ids_to_select))

# 准备用于标签的数据，并计算标签的固定 X 位置和对齐方式 (仅标签版本需要)
labeling_data <- NULL # 初始化为 NULL
label_x_limit_factor <- 1.15 # X轴扩展因子，为标签留空间
label_padding_factor <- 0.95 # 标签距离X轴扩展边缘的距离因子

# 确定轴范围
max_abs_x <- max(abs(combined_data$x_axis_val), na.rm = TRUE)
x_limit <- ceiling(max_abs_x / 5) * 5
x_limit <- max(x_limit, 10)
max_y <- max(combined_data$GeneRatioNumeric, na.rm = TRUE)
y_limit <- ceiling(max_y * 10) / 10
y_limit <- max(y_limit, 0.1)


if(nrow(highlighting_data) == 0) {
  warning("警告：根据提供的 ID 列表，没有找到任何可以高亮/标记的条目。请检查您的 ID 列表是否正确。")
} else {
  print(paste("找到", nrow(highlighting_data), "个条目用于高亮显示。"))
  # 计算标签位置 (仅当有需要高亮的条目时)
  labeling_data <- highlighting_data %>% # 从高亮数据开始准备标签数据
    arrange(Type, desc(GeneRatioNumeric)) %>%
    mutate(
      label_x_pos = ifelse(Type == "GO", x_limit * label_padding_factor, -x_limit * label_padding_factor),
      label_hjust = ifelse(Type == "GO", 1, 0)
    )
}

# --- 6. 定义绘图元素 ---

# 颜色
go_color <- "#E69F00"
kegg_color <- "#56B4E9"
plot_colors <- c("GO" = go_color, "KEGG" = kegg_color)
# 字体大小
base_font_size <- 12

# --- 创建基础绘图对象 (不含标签, 不含特定主题修改) ---
p_base <- ggplot() +
  # 层 1: 所有点
  geom_point(data = combined_data,
             aes(x = x_axis_val, y = GeneRatioNumeric, fill = Type, size = Count),
             shape = 21, alpha = 0.6, color = "grey80") +
  # 层 2: 高亮显示选中的点 (如果存在)
  { if(nrow(highlighting_data) > 0)
    geom_point(data = highlighting_data,
               aes(x = x_axis_val, y = GeneRatioNumeric, fill = Type, size = Count),
               shape = 21, color = "black", stroke = 0.8)
  } +
  # 标度
  scale_fill_manual(values = plot_colors, name = "Enrichment Type") +
  scale_size_continuous(name = "Gene Count", range = c(1.5, 7), breaks = c(0, 20, 40, 60)) +
  scale_x_continuous(
    # 为标签版本预留足够空间
    limits = c(-x_limit * label_x_limit_factor, x_limit * label_x_limit_factor),
    breaks = seq(-x_limit, x_limit, by = ifelse(x_limit > 50, 20, 10)),
    labels = abs,
    name = expression(paste("-log"[10], "(", italic("q"), "-value)"))
  ) +
  scale_y_continuous(
    limits = c(0, y_limit * 1.05),
    breaks = seq(0, y_limit, by = ifelse(y_limit > 0.5, 0.1, 0.05)),
    name = "Gene Ratio"
  ) +
  # 标题
  labs(
    title = paste0("GO/KEGG Enrichment Results for Common DEGs In MG ",dataset)
  ) +
  # 居中轴线
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
  # 图例顺序
  guides(fill = guide_legend(order = 1), size = guide_legend(order = 2))


# --- 定义特定的主题修改 (刻度线向内, 图例右上角) ---
theme_custom <- theme_classic(base_size = base_font_size) +
  theme(
    plot.title = element_text(hjust = 0.5, size = rel(1.1), face = "bold", margin = margin(b = 25)), # 增加标题下方边距
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(), # 隐藏默认轴线，因为我们手动绘制了中心线
    axis.title.x = element_text(size = rel(1.1), margin = margin(t = 10)),
    axis.title.y = element_text(size = rel(1.1), margin = margin(r = 10)),
    axis.text = element_text(color = "black", size = rel(0.9)),
    axis.ticks = element_line(color = "black", linewidth = 0.4),
    # --- Y轴刻度线向内 ---
    axis.ticks.length.y.left = unit(-0.15, "cm"), # 负值使刻度线向内绘制
    axis.ticks.length.x = unit(0.15, "cm"), # 保持X轴刻度线向外（或设为0）
    
    # --- 图例位置和样式 ---
    legend.position = c(0.98, 0.98), # 右上角坐标 (x, y from 0 to 1)
    legend.justification = c("right", "top"), # 对齐方式
    legend.box.just = "right",
    legend.margin = margin(t=0, r=0, b=0, l=0), # 减少图例边距
    legend.box.margin = margin(t=5, r=5, b=5, l=5), # 图例框边距
    legend.title = element_text(size = rel(0.9), face = "bold"),
    legend.text = element_text(size = rel(0.85)),
    legend.key.size = unit(0.5, 'cm'),
    legend.background = element_rect(fill="transparent", colour = NA), # 透明背景
    legend.box.background = element_rect(fill="transparent", colour = NA) # 透明框背景
  )


# --- 7. 生成无标签版本 ---
p_unlabeled <- p_base + theme_custom

# --- 8. 生成有标签版本 ---
p_labeled <- p_base +
  # 添加标签层 (如果 labeling_data 存在)
  { if (!is.null(labeling_data) && nrow(labeling_data) > 0)
    geom_text_repel(
      data = labeling_data,
      aes(x = label_x_pos, y = GeneRatioNumeric, label = Description, hjust = label_hjust),
      size = base_font_size * 0.28,
      color = "black", fontface = "plain",
      segment.color = NA, direction = "y", force = 0.5, force_pull = 0.1,
      box.padding = 0.15, point.padding = NA, max.overlaps = Inf,
      min.segment.length = Inf
    )
  } +
  # 应用主题
  theme_custom


# --- 9. 显示和保存图形 ---

# 显示无标签图
print(p_unlabeled)
# 保存无标签图
ggsave(paste0(dir_name,"GO_KEGG_common_NoLabels",dataset,".pdf"), plot = p_unlabeled, width = 7.5, height = 7, device = cairo_pdf)
ggsave(paste0(dir_name,"GO_KEGG_common_NoLabels",dataset,".png"), plot = p_unlabeled, width = 7.5, height = 7, dpi = 300, units = "in", bg = "white")

# 显示有标签图
print(p_labeled)
# 保存有标签图 (需要更宽的尺寸)
ggsave(paste0(dir_name,"GO_KEGG_common_SideLabels",dataset,".pdf"), plot = p_labeled, width = 9.5, height = 7, device = cairo_pdf) # 增加宽度
ggsave(paste0(dir_name,"GO_KEGG_common_SideLabels",dataset,".png"), plot = p_labeled, width = 9.5, height = 7, dpi = 300, units = "in", bg = "white") # 增加宽度