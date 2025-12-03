# ─────────────────────────────────────────────────────────────
# GLOBAL SETTINGS
# (请根据您的文件和路径设置)
# ─────────────────────────────────────────────────────────────
methods  <- "KEGG" # GO or KEGG

# 1. 数据集信息
cell1    <- 'LumHR'
cell2    <- 'LumHR'
cell3    <- 'LumHR'
dataset1 <- 'R-MG'
dataset2 <- 'S-MG'
dataset3 <- 'M-MG'

# 按照您原脚本的逻辑生成 Tag (用于读取文件)
dataset1_tag <- paste0(methods, "-result", "-", dataset1, "-", cell1)
dataset2_tag <- paste0(methods, "-result", "-", dataset2, "-", cell2)
dataset3_tag <- paste0(methods, "-result", "-", dataset3, "-", cell3)

# ★【修改点 1】更改文件名和标题的生成格式
# 格式: dataset1(cell1)-dataset2(cell2)-dataset3(cell3)
core_name  <- paste0(methods,'-',dataset1, "(", cell1, ")-", dataset2, "(", cell2, ")-", dataset3, "(", cell3, ")")

# 定义最终的标题和输出文件名（保持一致）
plot_title <- core_name
out_name   <- core_name  # 输出文件名不含后缀

top_n_lab  <- 20 # 只展示 Sum_Q 最高的 Top 20

# 路径设置
root_dir   <- "D:/Rfile/result/1202"

# 自动创建目录
dir.create(file.path(root_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# 定义输入文件 (保持原逻辑不变)
file_dataset1 <- file.path("D:/111/DEG-result", paste0(dataset1_tag, ".csv"))
file_dataset2 <- file.path("D:/111/DEG-result", paste0(dataset2_tag, ".csv"))
file_dataset3 <- file.path("D:/111/DEG-result", paste0(dataset3_tag, ".csv"))

# 定义输出文件路径
out_pdf       <- file.path(root_dir, "plots", paste0(out_name, ".pdf"))
out_csv       <- file.path(root_dir, "tables", paste0(out_name, ".csv"))

# ─────────────────────────────────────────────────────────────
# LOAD LIBRARIES
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(tools)
})

# ─────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# ─────────────────────────────────────────────────────────────
parse_ratio <- function(x) {
  ifelse(grepl("/", x),
         {parts <- str_split_fixed(x, "/", 2)
         as.numeric(parts[,1]) / as.numeric(parts[,2])},
         as.numeric(x))
}

load_go <- function(path, tag) {
  if (!file.exists(path)) stop("文件未找到: ", path)
  
  df <- read_csv(path, show_col_types = FALSE)
  
  sig_col <- intersect(c("p.adjust","qvalue","pvalue"), names(df))[1]
  if (is.na(sig_col)) stop("在文件中未找到 p.adjust/qvalue/pvalue: ", path)
  
  df %>%
    mutate(Description = str_remove(Description, " - Mus musculus \\(house mouse\\)")) %>%
    mutate(
      GeneRatio = parse_ratio(GeneRatio),
      q_log10   = -log10(.data[[sig_col]])
    ) %>%
    select(ID, Description, GeneRatio, q_log10) %>%
    rename_with(~paste0(tag, "_", .), c(GeneRatio, q_log10))
}

# ─────────────────────────────────────────────────────────────
# DATA PROCESSING
# ─────────────────────────────────────────────────────────────
message("Loading datasets...")
d1 <- load_go(file_dataset1, dataset1_tag)
d2 <- load_go(file_dataset2, dataset2_tag)
d3 <- load_go(file_dataset3, dataset3_tag)

message("Calculating 3-way Intersection...")

# 取三次交集
merged_shared <- d1 %>%
  inner_join(d2, by = c("ID", "Description")) %>%
  inner_join(d3, by = c("ID", "Description"))

message("共找到 ", nrow(merged_shared), " 个三者共有的通路。")

# 计算 Sum_Q 并排序
merged_shared <- merged_shared %>%
  mutate(Sum_Q = .data[[paste0(dataset1_tag, "_q_log10")]] + 
           .data[[paste0(dataset2_tag, "_q_log10")]] + 
           .data[[paste0(dataset3_tag, "_q_log10")]]) %>%
  arrange(desc(Sum_Q))

# 筛选 Top N 并保存表格
plot_data_wide <- merged_shared %>% slice_head(n = top_n_lab)
write_csv(plot_data_wide, out_csv)
message("交集表格已保存: ", out_csv)

# ─────────────────────────────────────────────────────────────
# RESHAPE & PLOT
# ─────────────────────────────────────────────────────────────
message("Preparing plot data...")

# 转长数据
plot_data_long <- plot_data_wide %>%
  select(-ID, -Sum_Q) %>%
  pivot_longer(
    cols = contains("_q_log10") | contains("_GeneRatio"),
    names_to = c("Group", ".value"),
    names_pattern = "^(.*)_(q_log10|GeneRatio)$" 
  )

# X轴显示优化：将长的 tag 映射回简洁的 dataset 名字
group_map <- setNames(c(dataset1, dataset2, dataset3), 
                      c(dataset1_tag, dataset2_tag, dataset3_tag))

plot_data_long <- plot_data_long %>%
  mutate(DisplayGroup = group_map[Group]) 

# 锁定 Description 顺序 (按 Sum_Q 倒序)
plot_data_long$Description <- factor(plot_data_long$Description, 
                                     levels = rev(plot_data_wide$Description))

# 锁定 X 轴顺序
plot_data_long$DisplayGroup <- factor(plot_data_long$DisplayGroup,
                                      levels = c(dataset1, dataset2, dataset3))

message("Generating Bubble Plot...")

p <- ggplot(plot_data_long, aes(x = DisplayGroup, y = Description)) +
  geom_point(aes(size = GeneRatio, color = q_log10)) +
  scale_color_gradient(low = "#3D5AA9", high = "#C83C3C", name = "-log10(q)") +
  scale_size_continuous(range = c(3, 8), name = "Gene Ratio") +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black", face = "bold"),
    axis.text.y = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid.major = element_line(linetype = "dashed", color = "grey90"),
    # ★【修改点 2】调小标题字号 (size = 10)
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10) 
  ) +
  labs(title = plot_title)

ggsave(out_pdf, p, width = 9, height = 1 + top_n_lab * 0.35, limitsize = FALSE)

message("✅ 图已保存：", out_pdf)