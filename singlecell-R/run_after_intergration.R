library(Seurat)
library(ggplot2)
library(dplyr)
library(limma)
library(cowplot)
library(wesanderson)
library(scales)
library(patchwork)
#########全局变量##########
dataset <- "MRS"

FP <- wes_palette("Zissou1", 5, type = "discrete")
palette <- wes_palette("Darjeeling1", length(sample_order), type = "continuous")

PC <- 20
res <- 0.4

# 初始化一个空的数据框来存储结果
expression_percentage_df <- data.frame()

execute <- FALSE  ##是否进行注释，第一次运行应设置为FALSE
##注释结果
NewCellType <- list("Luminal" = c(4,8,10),
                    "Basal" = c(2,5,6,12),
                    "NA1" = c(1,3,9,11),
                    "StemCells" = c(0,13,7)  )


#############################
updateCellType <- function(seurat_obj, cellTypeMapping, execute) {
  if (!execute) {
    return(seurat_obj)  # 如果execute为FALSE，则直接返回原Seurat对象
  }
  
  # 处理新的输入结构，将其转换为之前所需的结构
  newMapping <- unlist(lapply(names(cellTypeMapping), function(cellType) {
    setNames(rep(cellType, length(cellTypeMapping[[cellType]])), cellTypeMapping[[cellType]])
  }))
  newMapping <- newMapping[order(as.numeric(names(newMapping)))]
  # 确保cellTypeMapping向量包含了Seurat对象所有簇的映射
  all_clusters_present <- all(levels(Idents(seurat_obj)) %in% names(newMapping))
  if (!all_clusters_present) {
    stop("Not all clusters in the Seurat object are represented in the cellTypeMapping vector. Please include all clusters.")
  }
  
  # 更新Seurat对象的细胞类型
  seurat_obj$CellType <- newMapping[Idents(seurat_obj)]
  
  return(seurat_obj)  # 返回修改后的Seurat对象
}


################################
ea <- RunPCA(ea, features = VariableFeatures(object = ea))
ea <- FindNeighbors(ea, dims = 1:PC)
ea <- FindClusters(ea, resolution = res)
ea <- RunUMAP(ea, dims=1:PC)
ea <- updateCellType(ea, NewCellType, execute)
saveRDS(ea, paste0(dataset, "_CCAfinal.rds"))

p1_cca <- DimPlot(ea, reduction = "umap", group.by = "merge_dataset") +
  scale_colour_manual(values = palette) +
  ggtitle(paste0(dataset, " - CCA - By Dataset"))

p2_cca <- DimPlot(ea, reduction = "umap", group.by = "CellType") + ggtitle(paste0(dataset," - CCA - By CellType"))

p3_cca <- UMAPPlot(ea, label=T, label.size=8)+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  ggtitle(paste0(dataset," - CCA - By cluster"))

plotc1 <- p1_cca|p2_cca
plotc2 <- p1_cca|p2_cca|p3_cca
ggsave(paste0(dataset,"_CCA_celltype.png"), plot = plotc1, width = 20, height = 10)
ggsave(paste0(dataset,"_CCA_2.png"), plot = plotc2, width = 30, height = 10)
print("UMAP plots generated and saved.")

Idents(ea) <- ea$CellType


############CELLTYPE dotplot ######################
if (length(levels(ea@active.ident)) > 1) {
  markers <- FindAllMarkers(ea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 筛选每个cluster的top 5表达基因
  top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  features <- unique(as.character(top5_markers$gene))
  
  # 获取Zissou1配色方案
  palette <- wes_palette("Zissou1", 100, type = "continuous")
  
  # DotPlot 使用Zissou1配色方案
  p_dot <- DotPlot(ea, features = features) + RotatedAxis() +
    scale_colour_gradientn(colors = palette) # 应用配色方案
  
  # 保存DotPlot图像
  ggsave(paste0(dataset, "CCA_DotPlot.pdf"), plot = p_dot, width = 15, height = 8)
  ggsave(paste0(dataset, "CCA_celltype_DotPlot.tiff"), plot = p_dot, width = 15, height = 8, dpi = 300)
  
  
  # 继续你的DEG分析和保存逻辑
  write.csv(markers, paste0(dataset, "CCA_DEG.csv"))
}



