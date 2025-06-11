# =============================================================
# CellChat Visualization Pipeline (Optimized v3)
# =============================================================
# Author: ChatGPT
# Date: 2025‑05‑12
# -------------------------------------------------------------
# This script now generates, for every CellChat object:
#   1. Global interaction circle plots (count & weight)
#   2. Per‑cell‑type circle plots
#   3. Bubble plots of all ligand‑receptor interactions
#   4. Chord diagrams for selected pathways (with CSV exports)
#   5. River / Sankey plots of outgoing communication patterns
# All graphics are exported as high‑resolution PDF + PNG.
# =============================================================

# === 0. Environment === --------------------------------------------------------
library(CellChat)
library(patchwork)
library(ggplot2)
library(grid)        # annotate plots / titles
library(ggalluvial)  # river plots need this
library(NMF)
library(ggplotify)
# =============== USER SETTINGS ================================================
output_dir    <- "D:/111/CellChat"                   # change if needed
dataset_names <- c('R-MG'
                   #,'S-MG','M-MG','R-AG','S-AG','R-CG'
                   )
use_spefic_color <- TRUE
# Pathways for dedicated chord diagrams + CSV export
target_pathways <- c("RANKL", "WNT", "KIT", "DHT")
color_dict <- c(
  "ASC-rb"        = "#e31a1c",  # Strong Red
  "ASC-sg"        = "#a6cee3",  # Light Blue
  "Epi-Lgals7"    = "#808080",  # Gray
  "Basal"         = "#1f78b4",  # Strong Blue
  "CSC"           = "#fdbf6f",  # Light Orange/Peach
  "Epi-Fibro"     = "#cab2d6",  # Lavender
  "Epi-Pro"       = "#D3D3D3",  # LightGray
  "Lum-Immune"    = "#ffd700",  # Gold
  "Lum-Kit"       = "#008080",  # Teal
  "LumSEC-AG-Pip" = "#4B0082",  # Indigo
  "LumSEC-AG-t1"  = "#B8860B",  # DarkGoldenrod
  "LumSEC-AG-t2"  = "#fb9a99",  # Light Red/Pink
  "Epi-Krt7"      = "#b2df8a",  # Light Green
  "LumHR"         = "#ff7f00",  # Strong Orange
  "Lum-Basal"     = "#778899",  # LightSlateGray
  "Epi-Hmgcr"     = "#2F4F4F",  # DarkSlateGray
  "LumSEC-Lac"    = "#6a3d9a",  # Strong Purple
  "LumSEC-Lip"    = "#b15928",  # Brown
  "Epi-Lip"       = "#6B8E23",  # OliveDrab  (updated)
  "Epi-Lalba"     = "#87CEEB",  # SkyBlue
  "LumSEC-Mgp"    = "#00FF7F",  # SpringGreen
  "LumSEC-Vcam1"  = "#FF00FF",  # Magenta/Fuchsia
  "MaSC"          = "#33a02c",  # Strong Green
  "MaSC-Pro"      = "#FF69B4",  # HotPink
  "MaSC-t2-sg"    = "#FAEBD7",  # AntiqueWhite
  "LumSEC-Lip-CG" = "#D2B48C",  # Tan
  "Lum-Tm4sf4"    = "#00FFFF",  # Cyan/Aqua
  "Lum-Stat4"     = "#8FBC8F"   # DarkSeaGreen
)

# ==============================================================================
# 函数定义
filter_color_scheme <- function(color_dict, seurat_object) {
  existing_celltypes <- unique(cellchat@meta$newcelltype)
  
  filtered_colors <- color_dict[names(color_dict) %in% existing_celltypes]
  
  # 返回过滤后的配色向量
  return(filtered_colors)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: safe replay of current base‑graphics plot -------------------------
save_current_plot <- function(file_stub, width = 6, height = 6, res = 300) {
  rp <- recordPlot()                               # capture base‑graphics plot
  pdf(file.path(output_dir, paste0(file_stub, ".pdf")), width = width, height = height)
  replayPlot(rp)
  dev.off()
  png(file.path(output_dir, paste0(file_stub, ".png")), width = width * res,
      height = height * res, res = res)
  replayPlot(rp)
  dev.off()
}

# --- Helper: write communication table ----------------------------------------
write_comm_table <- function(cellchat, dataset, pathway) {
  tbl <- subsetCommunication(cellchat, signaling = pathway)
  write.csv(tbl, file.path(output_dir, sprintf("%s_%s_details.csv", dataset, pathway)),
            row.names = FALSE)
}

# === 1. Load CellChat objects ==================================================
cat("→ Loading CellChat objects…\n")
file_paths    <- file.path(output_dir, paste0(dataset_names, "_CellChat.rds"))
cellchat_list <- lapply(file_paths, function(fp) if (file.exists(fp)) readRDS(fp) else NULL)
names(cellchat_list) <- dataset_names

# === 2. Iterate over datasets ==================================================
for (dataset in dataset_names) {
  cellchat <- cellchat_list[[dataset]]
  if (is.null(cellchat)) {
    warning(dataset, " CellChat object is NULL – skipping…")
    next
  }
  if(use_spefic_color){
    filtered_colors <- filter_color_scheme(color_dict, cellchat) 
  }else{filtered_colors <- NULL}

  cat("\n========== Processing", dataset, "==========\n")
  grp_size <- as.numeric(table(cellchat@idents))
  mat      <- cellchat@net$weight
  edge_max <- max(mat)
  
  # ---------- 2.1 Global Circle Plots -----------------------------------------
  cat("  • Global circle plots…\n")
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count,  vertex.weight = grp_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Number of interactions",color.use = filtered_colors)
  netVisual_circle(cellchat@net$weight, vertex.weight = grp_size, weight.scale = TRUE,
                   label.edge = FALSE, title.name = "Interaction weights/strength",color.use = filtered_colors)
  save_current_plot(paste0(dataset, "_global_circle"))
  graphics.off()
  
  # ---------- 2.2 Per‑cell‑type Circle Plots ----------------------------------
  cat("  • Per‑cell‑type circle plots…\n")
  n_row <- nrow(mat)
  par(mfrow = c(ceiling(n_row / 4), 4), xpd = TRUE)
  for (i in seq_len(n_row)) {
    mat2 <- matrix(0, nrow = n_row, ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = grp_size, weight.scale = TRUE,
                     edge.weight.max = edge_max, title.name = rownames(mat)[i],color.use = filtered_colors)
  }
  save_current_plot(paste0(dataset, "_celltype_circles"), width = 8, height = 8)
  graphics.off()
  
  # ---------- 2.3 Bubble Plot --------------------------------------------------
  cat("  • Bubble plot…\n")
  bubble_plot <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL,
                                  remove.isolate = FALSE) +
    ggtitle(paste0(dataset, " – Cell‑Cell Communication Bubble Plot"))
  ggsave(file.path(output_dir, paste0(dataset, "_Communication_BubblePlot.png")),
         bubble_plot, width = 32, height = 24, dpi = 300)
  # optional PDF (comment out if unnecessary)
  ggsave(file.path(output_dir, paste0(dataset, "_Communication_BubblePlot.pdf")),
         bubble_plot, width = 32, height = 24)
  
  # ---------- 2.4 Pathway‑specific Chord Diagrams -----------------------------
  cat("  • Pathway‑specific chord diagrams…\n")
  available_pws <- intersect(target_pathways, cellchat@netP$pathways)
  for (pw in available_pws) {
    write_comm_table(cellchat, dataset, pw)
    par(mfrow = c(1,1))
    netVisual_aggregate(cellchat, signaling = pw, layout = "chord",color.use = filtered_colors)
    save_current_plot(paste0(dataset, "_", pw, "_chord"))
    graphics.off()
  }
}

cat("\nAll processing finished. Results saved to:", output_dir, "\n")

# ---------- 2.5 River Plot (Outgoing patterns) -----------------------------
cat("  • River / Sankey plot…\n")
dataset <- 'S-MG'
file_paths    <- file.path(output_dir, paste0(dataset, "_CellChat.rds"))
cellchat <- readRDS(file_paths)
if(use_spefic_color){
  filtered_colors <- filter_color_scheme(color_dict, cellchat) 
}else{filtered_colors <- NULL}
# Determine optimal k (selectK) then identify patterns
selectK(cellchat, pattern = "outgoing")
nPatterns <- 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing",color.use = filtered_colors)
save_current_plot(paste0(dataset, "_river"), width = 7, height = 5)
graphics.off()

