#############################################
## 1. Setup Environment and Load Packages
#############################################
library(CellChat)
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ggalluvial)
library(NMF)
library(ggplotify)
#############################################
## 2. Load Data and Check Cell Annotations
#############################################
# Replace "your_input_file.rds" with the path to your RDS file
dataset <- '23heart'
seuratObj <- readRDS(paste0("D:/111/combined_seurat_full_23_months_20k.rds"))
############################################
seuratObj <- NormalizeData(seuratObj)
assay <- seuratObj[["RNA"]]


# 1) 取 gene_name，缺失/空则用原 Ensembl 行名兜底（不改矩阵）
gene_name <- as.character(assay@meta.features$gene_name)
if (is.null(gene_name)) stop("RNA@meta.features$gene_name 不存在")
gene_name <- trimws(gene_name)

old_names <- rownames(assay)               # 原行名（Ensembl）
use_old   <- is.na(gene_name) | gene_name == ""
gene_name[use_old] <- old_names[use_old]   # 兜底：保留 Ensembl

# 2) Seurat 要求行名唯一：仅对“行名”做去重后缀，不改 gene_name 本身
new_names <- make.unique(gene_name, sep = "_dup")

# 建立映射：old -> new
map_old2new <- setNames(new_names, old_names)

# 3) 仅更新“行名”，不改矩阵数值
# counts / data / scale.data 的行名都换成 new_names（若存在）
if (!is.null(assay@counts) && nrow(assay@counts) > 0)
  rownames(assay@counts) <- new_names
if (!is.null(assay@data) && nrow(assay@data) > 0)
  rownames(assay@data) <- new_names
if (!is.null(assay@scale.data) && nrow(assay@scale.data) > 0) {
  # scale.data 的行名是特征名子集，这里用映射替换，保持一一对应
  rn <- rownames(assay@scale.data)
  rn_new <- ifelse(rn %in% names(map_old2new), map_old2new[rn], rn)
  rownames(assay@scale.data) <- rn_new
}

# 4) var.features 也改名（是字符向量）
if (length(assay@var.features)) {
  vf <- assay@var.features
  vf_new <- ifelse(vf %in% names(map_old2new), unname(map_old2new[vf]), vf)
  assay@var.features <- vf_new
}

# 5) 更新 meta.features：行名改为 new_names；保留 gene_id=原Ensembl，gene_name=你提供的名字
mf <- assay@meta.features
mf$gene_id   <- old_names           # 原始 Ensembl
mf$gene_name <- gene_name           # 原始想要的基因名（无 _dup 后缀）
rownames(mf) <- new_names
assay@meta.features <- mf

# 6) 回写对象
seuratObj[["RNA"]] <- assay
#############################################
## 3. Build the CellChat Object
#############################################
# Use the log-transformed expression matrix (stored in RNA@data by default)
data.input <- GetAssayData(seuratObj, assay = "RNA", slot = "data")
meta <- seuratObj@meta.data

# Build the CellChat object using "newcelltype" as the cell grouping information
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "Main_cell_type")

# Select the mouse version of the CellChat database (make sure CellChatDB.mouse is installed)
CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB

#############################################
## 4. Data Preprocessing and Communication Network Construction
#############################################
# Subset the data to only include genes of interest from the database
cellchat <- subsetData(cellchat)

# Identify overexpressed genes and ligand-receptor interactions in each cell group
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability between cells
cellchat <- computeCommunProb(cellchat)
# Filter out cell groups with too few cells (adjust the min.cells parameter if necessary)
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Compute the communication probability at the signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the ligand-receptor pairs to construct the overall communication network
cellchat <- aggregateNet(cellchat)

#############################################
## 5. Network Visualization
#############################################

# 5.1 Overall Cell-Cell Communication Network: Circular plot showing communication strength between cell groups
groupSize <- as.numeric(table(cellchat@idents))
overall_circle <- netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                                   weight.scale = TRUE, label.edge = FALSE,title.name = paste0(dataset, " - Overall Cell-Cell Communication Network"))
png(filename = paste0(dataset, "_Overall_Network_CirclePlot.png"), 
    width = 6, height = 6, units = "in", res = 300)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = TRUE, label.edge = FALSE)
dev.off()

# 5.2 Bubble Plot: Display all ligand-receptor interactions between cell groups
bubble_plot <- netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
bubble_plot <- bubble_plot + ggtitle(paste0(dataset, " - Cell-Cell Communication Bubble Plot"))
ggsave(paste0(dataset, "_Communication_BubblePlot.png"), plot = bubble_plot, width = 12, height = 8, dpi = 300)
print(bubble_plot)
p1 <- selectK(cellchat, pattern = "outgoing")
p1
nPatterns = 5
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
png(filename = paste0(dataset, "_Overall_Network_SankeyPlot.png"), 
    width = 6, height = 6, units = "in", res = 300)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#############################################
## 7. Additional Visualizations
#############################################


