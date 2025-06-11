library(Seurat)
library(harmony)
library(cowplot)
library(ggplot2)
library(dplyr)
library(limma)
library(wesanderson)
library(scales)
library(homologene)
###########定义全局变量###############
dataset <- "S-AG"
PCs <- 50
# 假设表达阈值为0.3
expression_threshold <-0.3
 # 设置随机数生成的种子以保证可重复性！！
set.seed(12345)
# 初始化一个空的数据框来存储结果
expression_percentage_df <- data.frame()
use_cluster <- FALSE # 是否基于cluster进行筛选
use_celltype <- TRUE #是否基于celltype筛选
sub_ident <- c() # 筛选的clusterID
celltype_list <- c("Basal cell","Luminal cell 1",'Luminal cell 2')# 筛选的CELLTYPE
cell_number <- 99999 # 每个cluster中想要随机选取的细胞数
execute <- FALSE ##是否进行注释，第一次运行应设置为FALSE
##注释结果
NewCellType <- list("Epith" = c(0,1),
"Luminal lact" = c(2,3),
"Luminal lact" = c(4,9,1)
)

#############################
random_subset_by_celltype <- function(seurat_obj, celltype_list, cell_number) {
  # 检查celltype_list是否为空
  if (length(celltype_list) == 0) {
    stop("celltype_list is empty. Please provide a list of cell types.")
  }
  
  # 检查celltype_list中的每个celltype是否存在于Seurat对象的CellType元数据中
  available_celltypes <- unique(seurat_obj$celltype)
  invalid_celltypes <- setdiff(celltype_list, available_celltypes)
  if (length(invalid_celltypes) > 0) {
    stop(paste("The following cell types are not present in the Seurat object:", paste(invalid_celltypes, collapse = ", "), ". Please provide valid cell types."))
  }
  
  # 根据CellType筛选Seurat对象
  seurat_subset <- subset(seurat_obj, subset = celltype %in% celltype_list)

  # 初始化一个空列表来存储每个celltype随机选出的细胞
  selected_cells_list <- list()

  # 对每个celltype进行循环处理
  for(cruuent_celltype in celltype_list) {
    # 获取当前celltype的所有细胞
    cells_in_celltype <- WhichCells(object = seurat_subset, expression = celltype == cruuent_celltype)

    # 检查当前celltype是否有足够的细胞供选择
    if (length(cells_in_celltype) == 0) {
      warning(paste("No cells found for cell type:", cruuent_celltype, ". Skipping this cell type."))
      next
    }
    
    # 如果当前celltype中的细胞数大于指定的细胞数，则随机选取指定数量的细胞
    if(length(cells_in_celltype) > cell_number) {
      selected_cells <- sample(cells_in_celltype, cell_number)
    } else {
      # 如果不足指定数量，则使用该celltype中的所有细胞
      selected_cells <- cells_in_celltype
    }

    # 将选出的细胞添加到列表中
    selected_cells_list[[cruuent_celltype]] <- selected_cells
  }

  # 将列表中的所有细胞合并为一个向量
  all_selected_cells <- unlist(selected_cells_list)

  # 根据选出的细胞创建一个新的Seurat对象
  final_subset <- subset(seurat_subset, cells = all_selected_cells)

  # 返回最终的Seurat对象
  return(final_subset)
}

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

  # 获取Seurat对象中的簇标识符
  cluster_ids <- Idents(seurat_obj)

  # 创建一个与Seurat对象细胞数量相同的字符向量，用于存储新的细胞类型
  new_celltype <- vector("character", length(cluster_ids))

  # 遍历所有簇标识符，更新new_celltype向量中的细胞类型
  for (cluster in unique(cluster_ids)) {
    new_celltype[cluster_ids == cluster] <- newMapping[as.character(cluster)]
  }

  # 更新Seurat对象的细胞类型
  seurat_obj$celltype <- new_celltype

  return(seurat_obj)  # 返回修改后的Seurat对象
}


#####################################
ea <- readRDS(paste0('../',dataset,".rds"))
ea$celltype <- ea$celltype_new ##annotation from h5ad
ea$oldcelltype <- ea$celltype_new ##annotation from h5ad
if(use_celltype){
   ea <- random_subset_by_celltype(ea, celltype_list, cell_number)
}

ea$batch_list <- ifelse(ea$rename == "S_GES_14_AG", 'batch1' , 'batch2')

obj <- ea
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$batch_list)
print("begin PCA")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj,npcs = PCs)
print("FINISH PCA")
obj <- IntegrateLayers(
  object = obj, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = "integrated.cca",
  verbose = FALSE
)
print("####################################### FINISH CCA ####################################")
obj <- FindNeighbors(obj, dims = 1:10,reduction = "integrated.cca")
obj <- FindClusters(obj, resolution = 0.2)
print("Running UMAP.")
obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:10)
print("UMAP completed.")
ea <- JoinLayers(obj)
ea <- updateCellType(ea, NewCellType, execute)
saveRDS(ea, paste0(dataset, "subset_for_monocle_CCA.rds"))
print("Seurat object saved.")
p1 <- DimPlot(ea, reduction = "umap", group.by = "rename") + 
  ggtitle(paste0(dataset, " - Subset - By Dataset"))

p2 <- DimPlot(ea, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1) + ggtitle(paste0(dataset," - Subset - By CellType - CCA"))

p3 <- UMAPPlot(ea, label=T, label.size=8)+
  theme(plot.title = element_text(size=20, face="bold"),
        axis.text=element_text(size=14),
        axis.title=element_text(size=16,face="bold"))+
  ggtitle(paste0(dataset," - Subset - By cluster"))
p4 <- DimPlot(ea, reduction = "umap", group.by = "time") +
  scale_color_viridis_d(option = 'D') +
  ggtitle("UMAP Plot colored by Time")
p5 <- DimPlot(ea, reduction = "umap", group.by = "oldcelltype", label = TRUE, pt.size = 1) + ggtitle(paste0(dataset," - Subset - By oldCellType"))
p6 <- DimPlot(ea, reduction = "umap", group.by = "gland", label = TRUE, pt.size = 1) + ggtitle(paste0(dataset," - Subset - By gland"))

plotc2 <- plotc2 <- (p1|p2|p3)/
	            (p4|p5|p6)

ggsave(paste0(dataset,"_Subset_celltype_CCA.png"), plot = plotc2, width = 30, height = 20)
ggsave(paste0(dataset,"_Subset_celltype.pdf"), plot = plotc2, width = 30, height = 20)

print("UMAP plots generated and saved.")


Idents(ea) <- ea$celltype
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

  ggsave(paste0(dataset, "celltype_DotPlot.tiff"), plot = p_dot, width = 15, height = 8, dpi = 300)
  # 继续你的DEG分析和保存逻辑
  write.csv(markers, paste0(dataset, "celltype_DEG.csv"))
}



Idents(ea) <- ea$seurat_clusters

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

  ggsave(paste0(dataset, "cluster_DotPlot.tiff"), plot = p_dot, width = 15, height = 8, dpi = 300)


  # 继续你的DEG分析和保存逻辑
  write.csv(markers, paste0(dataset, "cluster_DEG.csv"))
}






