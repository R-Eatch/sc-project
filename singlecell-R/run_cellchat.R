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
dataset <- 'R-MG'
seuratObj <- readRDS(paste0("D:/111/", dataset, "_for_cellchat.rds"))

#############################################
## 3. Build the CellChat Object
#############################################
# Use the log-transformed expression matrix (stored in RNA@data by default)
data.input <- GetAssayData(seuratObj, assay = "RNA", slot = "data")
meta <- seuratObj@meta.data

# Build the CellChat object using "newcelltype" as the cell grouping information
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "newcelltype")

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
ggsave(paste0(dataset, "_Communication_BubblePlot.png"), plot = bubble_plot, width = 32, height = 24, dpi = 300)
print(bubble_plot)
p1 <- selectK(cellchat, pattern = "outgoing")
p1
nPatterns = 6
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
png(filename = paste0(dataset, "_Overall_Network_SankeyPlot.png"), 
    width = 6, height = 6, units = "in", res = 300)
netAnalysis_river(cellchat, pattern = "outgoing")
dev.off()
#############################################
## 7. Additional Visualizations
#############################################


