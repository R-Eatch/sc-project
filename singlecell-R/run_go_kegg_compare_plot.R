# ─────────────────────────────────────────────────────────────
# GLOBAL SETTINGS  
# ─────────────────────────────────────────────────────────────
name <- '(Epi-Krt7)_GO_BP'   # ← 控制图名/输出路径等，可改为循环变量
dataset1_tag <- paste0("S-AG&R-AG", name)
dataset2_tag <- paste0("S-MG&S-AG", name)
plot_title   <- paste0("intersection_",name)  # 去除括号，更适合用作图标题
top_n_lab    <- 5

root_dir      <- "D:/Rfile/result/0620"
file_dataset1 <- file.path(root_dir, "tables", paste0(dataset1_tag, ".csv"))
file_dataset2 <- file.path(root_dir, "tables", paste0(dataset2_tag, ".csv"))
out_pdf       <- file.path(root_dir, "plots", paste0(plot_title, ".pdf"))
out_csv       <- file.path(root_dir, "tables", paste0(plot_title, ".csv"))
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(tidyr)
  library(ggplot2); library(ggrepel); library(tools)
})

# helper for -log10 parsing ------------------------------------------------------
parse_ratio <- function(x) {0
  ifelse(grepl("/", x),
         {parts <- str_split_fixed(x, "/", 2)
         as.numeric(parts[,1]) / as.numeric(parts[,2])},
         as.numeric(x))
}

load_go <- function(path, tag) {
  df <- read_csv(path, show_col_types = FALSE)
  sig_col <- intersect(c("p.adjust","qvalue","pvalue"), names(df))[1]
  if (is.na(sig_col)) stop("No p.adjust/qvalue/pvalue in ", path)
  df %>%
    mutate(Tissue    = tag,
           GeneRatio = parse_ratio(GeneRatio),
           q_log10   = -log10(.data[[sig_col]])) %>%
    select(ID, Description, Tissue, GeneRatio, q_log10)
}

d1 <- load_go(file_dataset1, dataset1_tag)
d2 <- load_go(file_dataset2, dataset2_tag)

merged <- full_join(
  d1 %>% rename(d1_q = q_log10, d1_ratio = GeneRatio),
  d2 %>% rename(d2_q = q_log10, d2_ratio = GeneRatio),
  by = c("ID","Description")
) %>%
  mutate(
    Category = case_when(
      !is.na(d1_q) &  is.na(d2_q) ~ dataset1_tag,
      is.na(d1_q)  & !is.na(d2_q) ~ dataset2_tag,
      TRUE                         ~ "Shared"
    ),
    PointRatio = coalesce(d1_ratio, d2_ratio),
    d1_q = replace_na(d1_q, 0),
    d2_q = replace_na(d2_q, 0),
    PointRatio = replace_na(PointRatio, 0)
  )

# ── save intersection CSV -------------------------------------------------------
if (out_csv == "")
  out_csv <- file.path(dirname(out_pdf),
                       paste0(dataset1_tag, "_VS_", dataset2_tag,
                              "_intersection.csv"))

merged %>%
  filter(Category == "Shared") %>%
  arrange(desc(d1_q + d2_q)) %>%
  select(ID, Description,
         d1_q, d2_q, d1_ratio, d2_ratio) %>%          # ① keep raw cols
  rename(                                             # ② give them new names
    !!dataset1_tag              := d1_q,
    !!dataset2_tag              := d2_q,
    !!paste0(dataset1_tag,
             "_ratio")          := d1_ratio,
    !!paste0(dataset2_tag,
             "_ratio")          := d2_ratio
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
if (plot_title == "")
  plot_title <- paste(dataset1_tag, "vs", dataset2_tag)

# ── PLOT ------------------------------------------------------------------------
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

ggsave(out_pdf, p, width = 10, height = 10)
message("✅ 图已保存：", out_pdf)
