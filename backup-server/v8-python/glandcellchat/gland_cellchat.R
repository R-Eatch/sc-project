#############################################
## 0. helper: choose DB & timestamp
#############################################
getSpeciesDB <- function(dataset) {
  # 根据 gland_pp.py 中的基因名格式 (如 Esr1, Epcam)，推测为小鼠数据
  # 如果是人类数据，请改为 CellChatDB.human
  CellChatDB.mouse 
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

# 输入目录需与 h5ad2seurat 脚本的输出目录一致
data_dir   <- "./data"
out_dir    <- "./fig_cellchat"
dir.create(out_dir, showWarnings = FALSE)

# 【修改1】数据集列表与 gland_pp.py 保持一致
dataset_list <- c('MG', 'SG', 'EG')

# 定义感兴趣的通路 (需确保这些通路在你的数据中存在)
focus_pathways <- c("RANKL", "WNT", "KIT", "DHT") 

# ++++ 全局随机抽样设置 (防止内存溢出) ++++
SUBSET_N_CELLS <- 20000
RANDOM_SEED    <- 1234

# ++++ 定义Seurat对象抽样函数 ++++
subsetSeuratObject <- function(seu_obj, n_target_cells, seed) {
  n_cells_original <- ncol(seu_obj)
  
  if (n_cells_original > n_target_cells) {
    message(ts(), "  Original cells: ", n_cells_original, ". Subsetting to ", n_target_cells, " random cells (seed=", seed, ")")
    set.seed(seed)
    sampled_cells <- sample(Cells(seu_obj), n_target_cells)
    seu_obj <- subset(seu_obj, cells = sampled_cells)
  } else {
    message(ts(), "  Using all ", n_cells_original, " cells (<= subset target)")
  }
  return(seu_obj)
}

#############################################
## 2. batch loop with try-catch
#############################################
for (ds in dataset_list) {
  cat("======", ds, "======\n")
  try({
    t0 <- ts(); message(t0, "  [", ds, "]  loading Seurat")
    
    # 【修改2】文件名匹配上一步转换脚本的输出 (_for_cellchat.rds)
    rds_file <- file.path(data_dir, paste0(ds, "_for_cellchat.rds"))
    
    if (!file.exists(rds_file)) {
      stop(paste("File not found:", rds_file))
    }
    
    seu <- readRDS(rds_file)
    
    # 执行抽样
   #seu <- subsetSeuratObject(seu, SUBSET_N_CELLS, RANDOM_SEED)
    seu <- NormalizeData(seu)
    # 提取数据
    data.input <- GetAssayData(seu, slot = "data")
    meta <- seu@meta.data
    
    # 【修改3】确保使用正确的细胞类型列名 (对应 gland_pp.py 中的 celltype)
    # 如果你的 metadata 里列名是 'leiden' 或其他，请在此处修改
    cell_type_column <- "celltype" 
    
    # 检查列是否存在，不存在则尝试回退到 'leiden'
    if (!cell_type_column %in% colnames(meta)) {
      if ("leiden" %in% colnames(meta)) {
        message("Warning: 'celltype' column not found, using 'leiden' instead.")
        cell_type_column <- "leiden"
      } else {
        stop(paste("Column", cell_type_column, "not found in metadata."))
      }
    }
    
    # 转换因子水平，确保顺序一致性（可选）
    meta[[cell_type_column]] <- as.factor(meta[[cell_type_column]])
    
    message(ts(), "  [", ds, "]  create CellChat using column: ", cell_type_column)
    chat <- createCellChat(data.input, meta = meta, group.by = cell_type_column)
    chat@DB <- getSpeciesDB(ds)
    
    # 标准处理流程
    chat <- subsetData(chat)
    chat <- identifyOverExpressedGenes(chat)
    chat <- identifyOverExpressedInteractions(chat)
    chat <- computeCommunProb(chat)
    
    # 过滤细胞数过少的交互
    chat <- filterCommunication(chat, min.cells = 10)
    
    chat <- computeCommunProbPathway(chat)
    chat <- aggregateNet(chat)
    
    ## 保存对象 (可选，解除注释以保存)
    # saveRDS(chat, file.path(out_dir, paste0(ds, "_CellChat.rds")))
    
    ## 绘图
    # 1. Circle Plot
    groupSize <- as.numeric(table(chat@idents))
    png(filename = file.path(out_dir, paste0(ds, "_Overall_Network_CirclePlot.png")), 
        width = 8, height = 8, units = "in", res = 300)
    par(mfrow = c(1,1))
    netVisual_circle(chat@net$weight, vertex.weight = groupSize, 
                     weight.scale = TRUE, label.edge = FALSE, 
                     title.name = paste0(ds, " Interaction Weights"))
    dev.off()
    pdf(
    file = file.path(out_dir, paste0(ds, "_Overall_Network_CirclePlot.pdf")),
    width = 6, height = 6, onefile = TRUE
    )
    netVisual_circle(
    chat@net$weight, vertex.weight = groupSize,
    weight.scale = TRUE, label.edge = FALSE,
    title.name = paste0(ds, " Interaction Weights")
    )
    dev.off()
    # 2. Bubble Plot
    bubble_plot <- netVisual_bubble(chat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
    bubble_plot <- bubble_plot + ggtitle(paste0(ds, " - Cell-Cell Communication Bubble Plot"))
    ggsave(file.path(out_dir, paste0(ds, "_Communication_BubblePlot.png")), 
           plot = bubble_plot, width = 36, height = 30, dpi = 300) # 调整了尺寸，防止过大
    ggsave(file.path(out_dir, paste0(ds, "_Communication_BubblePlot.pdf")),
           plot = bubble_plot, width = 36, height = 30) 
    # 3. Sankey/River Plot (Outgoing patterns)
    # 注意：如果 pattern 数多于实际识别出的模式，可能会报错，这里加个简单的 tryCatch 保护
    try({
      nPatterns = 5
      chat <- identifyCommunicationPatterns(chat, pattern = "outgoing", k = nPatterns)
      river_plot <- netAnalysis_river(chat, pattern = "outgoing")
      ggsave(file.path(out_dir, paste0(ds, "_Overall_Network_SankeyPlot.png")), 
             plot = river_plot, width = 10, height = 8, units = "in", dpi = 300)
      ggsave(file.path(out_dir, paste0(ds, "_Overall_Network_SankeyPlot.pdf")),
             plot = river_plot, width = 10, height = 8, units = "in")

    }, silent = TRUE)

    ## 导出表格
    # LR table
    lr_tab <- subsetCommunication(chat, slot.name = "net") %>% mutate(dataset = ds)
    write.csv(lr_tab, file.path(out_dir, paste0(ds, "_LR_table.csv")), row.names = FALSE)
    
    # Pathway table
    path_tab <- subsetCommunication(chat, slot.name = "netP") %>% mutate(dataset = ds)
    write.csv(path_tab, file.path(out_dir, paste0(ds, "_Pathway_table.csv")), row.names = FALSE)
    
    message(ts(), "  [", ds, "]  Finished successfully")
    
  }, silent = FALSE) -> err  # end try
  
  if (inherits(err, "try-error")) {
    logf <- file.path(out_dir, paste0(ds, "_error.log"))
    writeLines(c(ts(), err), logf)
    message(ts(), "  [", ds, "]  !!! Error logged to ", logf, "\n")
  }
}

#############################################
## 3. merge csv across datasets
#############################################
message(ts(), "  Merging results...")
lr_files   <- list.files(out_dir, "_LR_table.csv$", full.names = TRUE)
path_files <- list.files(out_dir, "_Pathway_table.csv$", full.names = TRUE)

if (length(lr_files) > 0) {
  do.call(rbind, lapply(lr_files, read.csv)) |>
    write.csv(file.path(out_dir, "All_LR_table.csv"), row.names = FALSE)
}
if (length(path_files) > 0) {
  do.call(rbind, lapply(path_files, read.csv)) |>
    write.csv(file.path(out_dir, "All_Pathway_table.csv"), row.names = FALSE)
}

message(ts(), "  === Batch CellChat Analysis Finished ===\n")
