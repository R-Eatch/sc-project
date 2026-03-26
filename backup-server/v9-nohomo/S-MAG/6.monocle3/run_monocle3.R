# Monocle3 pseudotime analysis with automated root selection via principal graph nodes

library(Seurat)
library(monocle3)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)

#### 定义全局变量 ###
dataset <- "S-MAG"
cores <- 8
####################

##### 读取数据并构建 CDS #####
ea <- readRDS(paste0(dataset, "_for_monocle.rds"))
data <- GetAssayData(ea, assay = 'RNA', slot = 'counts')
cell_metadata <- ea@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data), row.names = rownames(data))
cds <- new_cell_data_set(data, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
rm(data, cell_metadata, gene_annotation)

# 预处理与降维
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, cores = cores)
#plot_cells(cds, reduction_method = "UMAP", color_cells_by = "newcelltype")

# 使用 Seurat UMAP embedding
int.embed <- Embeddings(ea, reduction = "umap")
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed
#plot_cells(cds, reduction_method = "UMAP", color_cells_by = "newcelltype")

# 保存中间结果
save_monocle_objects(cds = cds, directory_path = paste0(dataset, '_cds_for_subset.rds'))

# 加载并筛选子集
cds_subset <- cds

# 聚类与学习主图
cds_subset <- cluster_cells(cds_subset)
cds_subset <- learn_graph(cds_subset, use_partition = FALSE)

# === 自动化选择根节点 (root_pr_nodes) ===
# 函数：根据 cell type 字段找出该类型细胞落在主图上最常见的 graph node
library(igraph)
get_root_pr_node <- function(cds_obj, celltype_col = 'newcelltype', target_type = 'StemCells') {
  # 找出目标类型的细胞索引
  cells <- colnames(cds_obj)[colData(cds_obj)[, celltype_col] == target_type]
  # 最近的 principal graph 节点映射
  closest_vert <- cds_obj@principal_graph_aux[['UMAP']]$pr_graph_cell_proj_closest_vertex
  closest_vert <- closest_vert[cells, , drop = FALSE]
  # 统计最频繁节点
  node_tab <- table(closest_vert)
  root_node <- names(node_tab)[which.max(node_tab)]
  return(root_node)
}

# 调用自动根节点选择
root_node <- get_root_pr_node(cds_subset, celltype_col = 'newcelltype', target_type = 'StemCells')
message('Selected root principal node: ', root_node)

# 使用 root_pr_nodes 进行 pseudotime 排序
cds_subset <- order_cells(cds_subset, root_pr_nodes = root_node)

# 保存最终 CDS
save_monocle_objects(cds = cds_subset, directory_path = paste0(dataset, '_subset_cds.rds'))

# === 绘图 ===
cds_subset <- load_monocle_objects(directory_path = paste0(dataset, '_subset_cds.rds'))

p1 <- plot_cells(cds_subset,
                 color_cells_by = 'pseudotime',
                 label_cell_groups = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 3)

p2 <- plot_cells(cds_subset,
                 color_cells_by = 'newcelltype',
                 label_groups_by_cluster = TRUE,
                 label_leaves = TRUE,
                 label_branch_points = TRUE,
                 graph_label_size = 3)

p1p2 <- wrap_plots(p1, p2)
ggsave(paste0(dataset, '_MONOCLE3.png'), p1p2, width = 15, height = 10)

ciliated_cds_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=cores)
pr_deg_ids <- subset(ciliated_cds_pr_test_res, q_value < 0.05)

Track_genes_sig <- pr_deg_ids %>% top_n(n=5, morans_I)%>% pull(gene_short_name) %>% as.character()


p3 <- plot_genes_in_pseudotime(cds_subset[Track_genes_sig,],
                               color_cells_by="newcelltype",
                               min_expr=0.5)
ggsave(paste0(dataset, "_Jitterplot.png"), plot = p3, width = 8, height = 6)

marker_test_res <- top_markers(cds_subset, group_cells_by="newcelltype", 
                               reference_cells=1000, cores=cores)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(5, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))


p4 <- plot_genes_by_group(cds,
                          top_specific_marker_ids,
                          group_cells_by="newcelltype",
                          ordering_type="maximal_on_diag",
                          max.size=5)

ggsave(paste0(dataset, "_dotplot1.png"), plot = p4, width = 8, height = 6)
