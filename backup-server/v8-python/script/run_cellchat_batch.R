#############################################
## 0. helper: choose DB & timestamp
#############################################
getSpeciesDB <- function(dataset) {
  # 简例：根据文件名前缀判定
  if (grepl("^M-|^R-", dataset)) {
    CellChatDB.mouse
  } else {
    CellChatDB.mouse  # 若无 marsupial DB，可先用 mouse 近似
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

data_dir   <- "D:/111/"
out_dir    <- "D:/111/CellChat/"
dir.create(out_dir, showWarnings = FALSE)

dataset_list <- c("R-MG", "R-AG", "S-MG", "S-AG", "M-MG","R-CG")
focus_pathways <- c("RANKL", "WNT", "KIT","DHT" )  # 通路名需与 DB 一致

#############################################
## 2. batch loop with try-catch
#############################################
for (ds in dataset_list) {
  cat("======", ds, "======\n")
  try({
    t0 <- ts(); message(t0, "  [", ds, "]  loading Seurat")
    seu <- readRDS(file.path(data_dir, paste0(ds, "_for_cellchat.rds")))
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
    saveRDS(chat, file.path(out_dir, paste0(ds, "_CellChat.rds")))
    message(ts(), "  [", ds, "]  CellChat object saved")
    
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
    chat.sub <- subsetCommunication(chat, signaling = focus_pathways)
    pdf(file.path(out_dir, paste0(ds, "_Bubble_Focus.pdf")),
        width = 8, height = 10)
    netVisual_bubble(chat.sub, angle.x = 45, remove.isolate = FALSE) +
      ggtitle(paste0(ds, "  focus: ", paste(focus_pathways, collapse = ", ")))
    dev.off()
    
    message(ts(), "  [", ds, "]  focus bubble saved\n")
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



