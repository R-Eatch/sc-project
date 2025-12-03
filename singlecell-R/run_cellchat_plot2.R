library(CellChat)
library(ggplot2)
dataset <- "S-MG"
cellchat <- readRDS(paste0("D:/111/CellChat/",dataset,"_CellChat.rds"))
focus_src <- "LumHR"

output_dir <- "D:/111/CellChat/"
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


## 1) 气泡图：从 Epi-krt7 指向所有 target（所有通路）
netVisual_bubble(
  cellchat,
  sources.use = focus_src,    # 也可给整数索引
  targets.use = NULL,         # NULL=所有 target；也可给向量
  signaling   = NULL,         # NULL=所有通路；也可 c("CDH","CLDN",...)
  remove.isolate = TRUE
)
save_current_plot(paste0(dataset,"_",focus_src,"_bubble_plot"))
# 文档: netVisual_bubble，点色=概率、点径=显著性/计数。:contentReference[oaicite:0]{index=0}

## 2) 聚合网络（按通路聚合后的“总强度”）+ 圈图/层级图
pathways.show <- cellchat@netP$pathways  # 或自己挑选若干通路
netVisual_aggregate(
  cellchat,
  signaling = pathways.show,  # 可给 1~数条通路
  layout = "circle",          # "circle"/"chord"/"hierarchy"/"spatial"
  sources.use = focus_src,    # 仅展示从 Epi-krt7 发出的边
  remove.isolate = TRUE
)
save_current_plot(paste0(dataset,"_",focus_src,"_circle_plot"))
# 文档: netVisual_aggregate 支持 sources.use / targets.use、layout 等。:contentReference[oaicite:1]{index=1}

## 3) 单条通路的“细粒度”圈图（建议挑 Top 通路逐一画）
netVisual_chord_cell(
  cellchat,
  signaling   = "CDH",        # 举例单通路
  sources.use = focus_src,    # 只看 Epi-krt7 发出的边
  lab.cex = 0.7               # 字体大小
)
# 文档: netVisual_chord_cell / netVisual_chord_gene。:contentReference[oaicite:2]{index=2}

## 4) outgoing 角色热图（行=通路，列=细胞群；值为相对强度）
netAnalysis_signalingRole_heatmap(
  cellchat,
  pattern = "outgoing"        # 或 "incoming"
)
# 该热图是“行归一化”的相对强度；你可在图上重点关注 Epi-krt7 这列。:contentReference[oaicite:3]{index=3}

## 5) 看某条通路内部由哪些 LR 对贡献最大
netAnalysis_contribution(cellchat, signaling = "CDH")
# 文档: netAnalysis_contribution。:contentReference[oaicite:4]{index=4}
