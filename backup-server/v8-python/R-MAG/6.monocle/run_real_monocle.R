library(Seurat)
library(monocle)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)
library(argparse)
library(viridis)
###################
##########################
####GlobalEnv example###
dataset <- "R-MAG"
cores <- 8
s.genes <- c("Esr1","Top2a","Tcf4")
brunch_gene <- s.genes
subset_celltype <- c("LumHR","StemCells","LumSEC")
root_state <- 1 #for ordercells
subset_clusters <- c()
branch_point <- 1 ##for BEAM analysis
set.seed(12345)

####################


#####定义函数#####

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
##################


ea <- readRDS(paste0(dataset,"_for_monocle.rds"))
ea$celltype=ea$newcelltype
subset_celltype <- unique(ea$newcelltype)
#ea <- random_subset_by_celltype(ea,subset_celltype,999999)
#ea <- subset(ea1, subset = seurat_clusters == subset_clusters)##基于cluster进行筛选
# Seurat to CDS
data <- GetAssayData(ea, assay = 'RNA', slot = 'counts') 
cell_metadata <- ea@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(data)) 
rownames(gene_annotation) <- rownames(data) 
pd <- new('AnnotatedDataFrame', data = cell_metadata)
fd <- new('AnnotatedDataFrame', data = gene_annotation)
# trajectory analysis
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds) 
cds <- estimateDispersions(cds)

express_genes <-rownames(ea[["RNA"]]@meta.features)[ea[["RNA"]]@meta.features$highly.variable == TRUE]
diff_test_res <- differentialGeneTest(cds[express_genes,],
                                      fullModelFormulaStr = "~sample")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
p11 <- plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = root_state)####定义根的起点
p1 <- plot_cell_trajectory(cds, color_by = "newcelltype", size = 10) 
p2 <- plot_cell_trajectory(cds,color_by="Pseudotime", size=10,show_backbone=TRUE)
p3 <- plot_cell_trajectory(cds, color_by = "stage") + scale_color_viridis_d(option = 'D')
p3.5 <- plot_cell_trajectory(cds, color_by = "State")
p3.6 <- plot_cell_trajectory(cds, color_by = "leiden",cell_name_size=10)
p4 <- plot_genes_jitter(cds[s.genes,], grouping = "newcelltype", color_by = "newcelltype")+
  ggtitle(paste0(dataset, " Gene Jitterplot")) +
  theme(
    legend.title = element_text(size = 20),  # 调整图例标题的大小
    legend.text = element_text(size = 20),   # 调整图例文本的大小
    plot.title = element_text(size = 20, hjust = 0.5)  # 调整标题的大小并居中
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5)))  # 调整图例中点的大小
p5 <- plot_genes_violin(cds[s.genes,], grouping = "newcelltype", color_by = "newcelltype")
p6 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "newcelltype")
plotc1 <- (p1 | p2) /
          (p3 | p3.5)
plotc2 <- p4|p5|p6
ggsave(paste0(dataset,"_pseudotime.png"), plot = plotc1, width = 20, height = 20)
ggsave(paste0(dataset,"_pseudotime.pdf"), plot = plotc1, width = 20, height = 20)
ggsave(paste0(dataset,"_Genes_Jitterplot.png"), plot = plotc2, width = 16, height = 8)
ggsave(paste0(dataset,"_cluster.png"), plot = p3.6, width = 10, height = 10)
ggsave(paste0(dataset,"_ordergenes.png"), plot = p11, width = 10, height = 10)
p3.7 <- plot_cell_trajectory(cds, color_by = "sample",cell_name_size=10)
ggsave(paste0(dataset,"_dataset.png"), plot = p3.7, width = 10, height = 10)

if(TRUE){

##BEAM analysis
BEAM_res <- BEAM(cds[ordering_genes,], branch_point = branch_point, cores = cores,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_GENE <- subset(BEAM_res, qval<0.01) %>% pull(gene_short_name) %>% as.character()
write.csv(BEAM_res, paste0(dataset,"Branch_all.csv"), row.names = F)

top100 <- top_n(BEAM_res, n=100, desc(qval)) %>% pull(gene_short_name) %>% as.character()

p8 <- plot_genes_branched_heatmap(cds[top100,],
                            branch_point = branch_point,
                            num_clusters = 3,
                            cores = cores,return_heatmap = T,
                            use_gene_short_name = T,
                            show_rownames = T)
ggsave(paste0(dataset,"_Branch_heatmap.png"), p8$ph_res, width = 6.5, height = 10)
df1 <- p8$annotation_row
df3 <- data.frame("gene" = rownames(df1), "cluster" = df1$Cluster)
write.csv(df3,paste0(dataset,"_Branch_cluster.csv"),row.names = F)
saveRDS(p8,"p8.rds")

intrest_gene <- row.names(subset(fData(cds),
                               gene_short_name %in% brunch_gene))
p9 <- plot_genes_branched_pseudotime(cds[intrest_gene,],
                               branch_point = branch_point,
                               color_by = "newcelltype",
                               ncol = 1)+
  ggtitle(paste0(dataset, " Branch Jitterplot"))
  guides(colour = guide_legend(override.aes = list(size = 5)))  # 调整图例中点的大小
ggsave(paste0(dataset,"_Branch_Jitterplot.png"), p9, width = 6.5, height = 10)
}

if(run <- TRUE){
  diff_test_res1 <- differentialGeneTest(cds[ordering_genes,],
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
  write.csv(diff_test_res1, paste0(dataset,"time_test.csv"), row.names = F)                                      
  top100 <- top_n(diff_test_res1, n=500, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  p <- plot_pseudotime_heatmap(cds[top100,],
                               num_clusters = 5,
                               cores = cores,return_heatmap=T,
                               show_rownames = F)
  
  clusters <- cutree(p$tree_row, k = 5) 
  clustering <- data.frame(clusters) 
  clustering[,1] <- as.character(clustering[,1]) 
  colnames(clustering) <- "Gene_Clusters" 
  table(clustering)
  write.csv(clustering, paste0(dataset,"_gene_Time_all.csv"), row.names = T)
  top100 <- top_n(diff_test_res1, n=500, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  p7 <- plot_pseudotime_heatmap(cds[top100,],
                                num_clusters = 5,
                                cores = cores,return_heatmap=T,
                                show_rownames = F)
  ggsave(paste0(dataset,"_Gene_heatmap.png"), p7, width = 6.5, height = 10)
ggsave(paste0(dataset,"_Gene_heatmap.pdf"), p7, width = 6.5, height = 10)

}
saveRDS(cds,paste0(dataset,"_monocle.rds"))


