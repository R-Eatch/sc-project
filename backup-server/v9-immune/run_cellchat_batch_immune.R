#############################################
## 0. helper: choose DB & timestamp
#############################################
getSpeciesDB <- function(dataset) {
  # 简例：根据文件名前缀判定
  if (grepl("^M-|^R-", dataset)) {
    CellChatDB.mouse
  } else {
    CellChatDB.mouse # 若无 marsupial DB，可先用 mouse 近似
  }
}

ts <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

#############################################
## 1. global settings
#############################################
library(CellChat)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggalluvial)

data_dir   <- "./forcellchat"
out_dir    <- "./result"
dir.create(out_dir, showWarnings = FALSE)

dataset_list <- c("M-MG",
                  "R-MG", "S-MG","S-AG", "R-AG","R-CG"
                  )
focus_pathways <- c("RANKL", "WNT", "KIT","DHT" ) # 通路名需与 DB 一致

# ++++ 新增：全局随机抽样设置 ++++
SUBSET_N_CELLS <- 20000
RANDOM_SEED    <- 1234

# ++++ 新增：定义Seurat对象抽样函数 ++++
#' @param seu_obj Seurat对象
#' @param n_target_cells 目标细胞数
#' @param seed 随机数种子
#' @return 经过抽样（如果需要）的Seurat对象
subsetSeuratObject <- function(seu_obj, n_target_cells, seed) {
  n_cells_original <- ncol(seu_obj)
  
  # 只有在当前细胞数 > 目标细胞数时才执行抽样
  if (n_cells_original > n_target_cells) {
    message(ts(), "  Original cells: ", n_cells_original, ". Subsetting to ", n_target_cells, " random cells (seed=", seed, ")")
    # 设定随机数种子以保证结果可重复
    set.seed(seed)
    # 随机抽取细胞名称
    sampled_cells <- sample(Cells(seu_obj), n_target_cells)
    # subset对象
    seu_obj <- subset(seu_obj, cells = sampled_cells)
  } else {
    # 如果细胞数本来就少于或等于目标，则不进行抽样
    message(ts(), "  Using all ", n_cells_original, " cells (<= subset target of ", n_target_cells, ")")
  }
  
  # 返回处理后的Seurat对象
  return(seu_obj)
}
# +++++++++++++++++++++++++++++++++++++


#############################################
## 2. batch loop with try-catch
#############################################
for (ds in dataset_list) {
  cat("======", ds, "======\n")
  try({
    t0 <- ts(); message(t0, "  [", ds, "]  loading Seurat")
    seu <- readRDS(file.path(data_dir, paste0(ds, "_imm_for_cellchat.rds")))
    
    # ++++ 修改：调用函数进行抽样 ++++
    # 使用刚刚定义的函数来处理seu对象
    seu <- subsetSeuratObject(seu, SUBSET_N_CELLS, RANDOM_SEED)
    # ++++++++++++++++++++++++++++++++
    
    data.input <- GetAssayData(seu, slot = "data")
    meta <- seu@meta.data
    
    message(ts(), "  [", ds, "]  create CellChat")
    chat <- createCellChat(data.input, meta = meta, group.by = "newcelltype")
    chat@DB <- getSpeciesDB(ds)
    
    chat <- subsetData(chat)
    chat <- identifyOverExpressedGenes(chat)
    chat <- identifyOverExpressedInteractions(chat)
    chat <- computeCommunProb(chat)
    chat <- filterCommunication(chat, min.cells = 10)
    chat <- computeCommunProbPathway(chat)
    chat <- aggregateNet(chat)
    
    ## save object
    #saveRDS(chat, file.path(out_dir, paste0(ds, "_CellChat.rds")))
    #message(ts(), "  [", ds, "]  CellChat object saved")
    
    ## plot
    groupSize <- as.numeric(table(chat@idents))
    overall_circle <- netVisual_circle(chat@net$weight, vertex.weight = groupSize,
                                       weight.scale = TRUE, label.edge = FALSE,title.name = paste0(ds, " - Overall Cell-Cell Communication Network"))
    
    
    png(filename = file.path(out_dir, paste0(ds, "_Overall_Network_CirclePlot.png")), 
        width = 6, height = 6, units = "in", res = 300)
    netVisual_circle(chat@net$weight, vertex.weight = groupSize, 
                     weight.scale = TRUE, label.edge = FALSE)
    dev.off()
    pdf(
    file = file.path(out_dir, paste0(ds, "_Overall_Network_CirclePlot.pdf")),
    width = 6, height = 6, onefile = TRUE
    )
    netVisual_circle(
    chat@net$weight, vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE,
    title.name = circle_title
    )
    dev.off()
    # 5.2 Bubble Plot: Display all ligand-receptor interactions between cell groups
    bubble_plot <- netVisual_bubble(chat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
    bubble_plot <- bubble_plot + ggtitle(paste0(ds, " - Cell-Cell Communication Bubble Plot"))
    
    # [!!!! 已修复 2 !!!!]
    ggsave(file.path(out_dir, paste0(ds, "_Communication_BubblePlot.png")), plot = bubble_plot, width = 45, height = 30, dpi = 300)
    print(bubble_plot)
    #p1 <- selectK(chat, pattern = "outgoing")
    #p1
    nPatterns = 5
    chat <- identifyCommunicationPatterns(chat, pattern = "outgoing", k = nPatterns)
    netAnalysis_river(chat, pattern = "outgoing")
    river_plot <- netAnalysis_river(chat, pattern = "outgoing")
    ggsave(file.path(out_dir, paste0(ds, "_Overall_Network_SankeyPlot.png")),plot = river_plot,
        width = 6, height = 10, units = "in", dpi = 300) 
    # [!!!! 已修复 3 !!!!]
    #png(filename = file.path(out_dir, paste0(ds, "_Overall_Network_SankeyPlot.png")), 
    #    width = 6, height = 10, units = "in", res = 300)
    #netAnalysis_river(chat, pattern = "outgoing")
    dev.off()
    
    ## export LR table
    lr_tab <- subsetCommunication(chat, slot.name = "net") %>%
      mutate(dataset = ds)
    write.csv(lr_tab,
              file.path(out_dir, paste0(ds, "_LR_table.csv")),
              row.names = FALSE)
    
    ## export pathway table
    path_tab <- subsetCommunication(chat, slot.name = "netP") %>%
      mutate(dataset = ds)
    write.csv(path_tab,
              file.path(out_dir, paste0(ds, "_Pathway_table.csv")),
              row.names = FALSE)
    
    message(ts(), "  [", ds, "]  CSV exported")
    
    ## bubble plot for focus pathways
    # chat.sub <- subsetCommunication(chat, signaling = focus_pathways)
    # pdf(file.path(out_dir, paste0(ds, "_Bubble_Focus.pdf")),
    #     width = 8, height = 10)
    # netVisual_bubble(chat.sub, angle.x = 45, remove.isolate = FALSE) +
    #   ggtitle(paste0(ds, "  focus: ", paste(focus_pathways, collapse = ", ")))
    # dev.off()
    
    # message(ts(), "  [", ds, "]  focus bubble saved\n")
  }, silent = FALSE) -> err  # end try
  
  if (inherits(err, "try-error")) {
    logf <- file.path(out_dir, paste0(ds, "_error.log"))
    writeLines(c(ts(), err), logf)
    message(ts(), "  [", ds, "]  !!! error logged to ", logf, "\n")
  }
}

#############################################
## 3. merge csv across datasets (optional)
#############################################
lr_files   <- list.files(out_dir, "_LR_table.csv$", full.names = TRUE)
path_files <- list.files(out_dir, "_Pathway_table.csv$", full.names = TRUE)

if (length(lr_files)) {
  do.call(rbind, lapply(lr_files, read.csv)) |>
    write.csv(file.path(out_dir, "All_LR_table.csv"), row.names = FALSE)
}
if (length(path_files)) {
  do.call(rbind, lapply(path_files, read.csv)) |>
    write.csv(file.path(out_dir, "All_Pathway_table.csv"), row.names = FALSE)
}

message(ts(), "  === batch finished ===\n")
