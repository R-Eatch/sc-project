# ─────────────────────────────────────────────────────────────
# GLOBAL SETTINGS  
# (请根据您的文件和路径设置)
# ─────────────────────────────────────────────────────────────
methods <- "GO"#GO or KEGG
cell1 <- 'LumHR'
cell2 <- 'LumSEC-AG-t2'
dataset1 <- 'R-MG'
dataset2 <- 'R-AG'


name <- paste0(methods,"_",dataset1,"-",cell1,"_",dataset2,"-",cell2)
dataset1_tag <- paste0(methods,"-result","-",dataset1,"-",cell1)
dataset2_tag <- paste0(methods,"-result","-",dataset2,"-",cell2)

# 自动生成更简洁的标题
plot_title   <- paste0("Intersection_", name) 

top_n_lab    <- 5 # 每组显示 top 5 个标签

# 【重要】请修改为您的实际路径
root_dir     <- "D:/Rfile/result/1106"
# 自动创建目录（如果不存在）
dir.create(file.path(root_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(root_dir, "plots"), recursive = TRUE, showWarnings = FALSE)

# 定义输入输出文件
file_dataset1 <- file.path("D:/111/DEG-result", paste0(dataset1_tag, ".csv"))
file_dataset2 <- file.path("D:/111/DEG-result", paste0(dataset2_tag, ".csv"))
out_pdf       <- file.path(root_dir, "plots", paste0(plot_title, ".pdf"))
out_csv       <- file.path(root_dir, "tables", paste0(plot_title, ".csv"))

# ─────────────────────────────────────────────────────────────
# LOAD LIBRARIES
# ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(tools)
})

# ─────────────────────────────────────────────────────────────
# HELPER FUNCTIONS
# (您提供的原始函数)
# ─────────────────────────────────────────────────────────────

# helper for -log10 parsing ------------------------------------------------------
parse_ratio <- function(x) {
  ifelse(grepl("/", x),
         {parts <- str_split_fixed(x, "/", 2)
         as.numeric(parts[,1]) / as.numeric(parts[,2])},
         as.numeric(x))
}

# ★ 修改后的 load_go 函数 ★
# ------------------------------------------------------------------------------
load_go <- function(path, tag) {
  # 检查文件是否存在
  if (!file.exists(path)) {
    stop("文件未找到: ", path)
  }
  
  df <- read_csv(path, show_col_types = FALSE)
  
  # 寻找显著性列 (p.adjust, qvalue, pvalue)
  sig_col <- intersect(c("p.adjust","qvalue","pvalue"), names(df))[1]
  if (is.na(sig_col)) stop("在文件中未找到 p.adjust/qvalue/pvalue: ", path)
  
  df %>%
    # ★ 新增步骤：清理KEGG描述 ★
    #   自动删除 " - Mus musculus (house mouse)" 字符串
    mutate(Description = str_remove(Description, " - Mus musculus \\(house mouse\\)")) %>%
    
    # ↓ 以下是您原来的处理步骤
    mutate(Tissue    = tag,
           GeneRatio = parse_ratio(GeneRatio),
           q_log10   = -log10(.data[[sig_col]])) %>%
    select(ID, Description, Tissue, GeneRatio, q_log10)
}

# ─────────────────────────────────────────────────────────────
# DATA PROCESSING
# (您提供的原始逻辑)
# ─────────────────────────────────────────────────────────────
message("Loading dataset 1: ", file_dataset1)
d1 <- load_go(file_dataset1, dataset1_tag)

message("Loading dataset 2: ", file_dataset2)
d2 <- load_go(file_dataset2, dataset2_tag)

message("Merging datasets...")
merged <- full_join(
  d1 %>% rename(d1_q = q_log10, d1_ratio = GeneRatio),
  d2 %>% rename(d2_q = q_log10, d2_ratio = GeneRatio),
  by = c("ID","Description")
) %>%
  mutate(
    # 自动分类：交集或独有
    Category = case_when(
      !is.na(d1_q) &  is.na(d2_q) ~ dataset1_tag,
      is.na(d1_q) & !is.na(d2_q) ~ dataset2_tag,
      TRUE                       ~ "Shared"
    ),
    # 合并 Ratio 和 q 值，用于绘图
    PointRatio = coalesce(d1_ratio, d2_ratio),
    d1_q = replace_na(d1_q, 0),
    d2_q = replace_na(d2_q, 0),
    PointRatio = replace_na(PointRatio, 0)
  )

# ── save intersection CSV -------------------------------------------------------
# (您提供的原始逻辑，用于输出交集表格)
message("Saving intersection CSV...")
merged %>%
  filter(Category == "Shared") %>%
  arrange(desc(d1_q + d2_q)) %>%
  select(ID, Description,
         d1_q, d2_q, d1_ratio, d2_ratio) %>%
  rename(
    !!dataset1_tag := d1_q,
    !!dataset2_tag := d2_q,
    !!paste0(dataset1_tag, "_ratio") := d1_ratio,
    !!paste0(dataset2_tag, "_ratio") := d2_ratio
  ) %>%
  write_csv(out_csv)

message("交集通路已保存：", out_csv)

# ── label subset ----------------------------------------------------------------
label_df <- merged %>%
  group_by(Category) %>%
  arrange(desc(d1_q + d2_q)) %>%
  slice_head(n = top_n_lab) %>%
  ungroup()

# ── colours & captions ----------------------------------------------------------
cols <- c(setNames("#3D5AA9", dataset1_tag),
          setNames("#C83C3C", dataset2_tag),
          Shared = "#7DA34D")

x_lab <- paste0("-log10(q)  [", dataset1_tag, "]")
y_lab <- paste0("-log10(q)  [", dataset2_tag, "]")

# ─────────────────────────────────────────────────────────────
# PLOT
# (您提供的原始绘图代码)
# ─────────────────────────────────────────────────────────────
message("Generating plot...")
p <- ggplot(merged, aes(d1_q, d2_q)) +
  geom_point(aes(color = Category, size = PointRatio),
             position = position_jitter(width=.08,height=.08),
             alpha=.85, stroke=.25) +
  geom_text_repel(data=label_df,
                  aes(label=Description, color=Category),
                  force=2, max.iter=20000,
                  box.padding=.5, point.padding=.3,
                  segment.size=.2, size=3, show.legend=FALSE) +
  scale_size_continuous(name = "Gene ratio", range=c(3,9),
                        limits=c(0,0.25),
                        breaks=c(.05,.10,.15,.20),
                        labels=scales::number_format(accuracy=.01)) +
  scale_color_manual(values = cols, breaks = c(dataset1_tag,
                                               dataset2_tag,"Shared")) +
  scale_x_continuous(expand=expansion(mult=.05), name=x_lab) +
  scale_y_continuous(expand=expansion(mult=.05), name=y_lab) +
  coord_cartesian() +
  labs(title = plot_title, color = "Dataset") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(colour="grey75", linewidth=.3),
        panel.grid.minor = element_blank(),
        legend.position  = "right",
        legend.box       = "vertical")

ggsave(out_pdf, p, width = 10, height = 10, bg = "white")
message("✅ 图已保存：", out_pdf)