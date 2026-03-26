library(Seurat)
library(monocle3)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)
####定义全局变量###
dataset <- "S-MAG"
cores <- 4
####################


#####定义函数#####

get_earliest_principal_node <- function(cds){
  # 不再需要筛选特定时间区间的细胞，所以我们取所有的细胞ID
  cell_ids <- colnames(cds)
  
  # 计算每个细胞最接近的主图谱节点
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[cell_ids, ])
  
  # 确定所有细胞中最常被映射到的主图谱节点
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex))))]
  
  return(root_pr_nodes)
}

##################


######运行如下代码######
ea <- readRDS(paste0("../6.monocle/",dataset,"_for_monocle.rds"))

data <- GetAssayData(ea, assay = 'RNA', slot = 'counts') 
cell_metadata <- ea@meta.data 
gene_annotation <- data.frame(gene_short_name = rownames(data)) 
rownames(gene_annotation) <- rownames(data) 
cds <- new_cell_data_set(data,cell_metadata = cell_metadata,gene_metadata = gene_annotation)

rm(data,cell_metadata,gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds,cores = cores)
plot_cells(cds, reduction_method="UMAP", color_cells_by="newcelltype")
int.embed <- Embeddings(ea, reduction = "umap")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="newcelltype")
save_monocle_objects(cds=cds, directory_path=paste0(dataset,'_cds_for_subset.rds'))
####marker test#########
marker_test_res <- top_markers(cds, group_cells_by="CellType", 
                               reference_cells=1000, cores=cores)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


p4 <- plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="CellType",
                    ordering_type="maximal_on_diag",
                    max.size=5)

ggsave(paste0(dataset, "_dotplot1.png"), plot = p4, width = 8, height = 6)

#####trajectory##############


cds <- cluster_cells(cds)
cds <- learn_graph(cds)

cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
p3 <- plot_cells(cds,
             color_cells_by = "pseudotime",
               label_cell_groups=TRUE,
             label_leaves=TRUE,
             label_branch_points=TRUE,
             graph_label_size=1.5)
ggsave(paste0(dataset, "_pse.png"), plot = p3, width = 8, height = 6)
