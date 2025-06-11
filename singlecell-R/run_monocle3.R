library(Seurat)
library(monocle3)
library(patchwork)
library(dplyr)
library(plotly)
library(ggplot2)
####定义全局变量###
dataset <- "M-MG"
cores <- 1 
####################


#####定义函数#####


##################


######运行如下代码######
#ea <- readRDS(paste0("D:/111/",dataset,"_for_monocle.rds"))
if(!file.exists(paste0(dataset,'_cds_for_subset.rds'))){
  ea <- readRDS(paste0(dataset,"_for_monocle.rds"))
  
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
}

#############电脑上操作#############
if(FALSE){
  cds <- load_monocle_objects(directory_path=paste0('D:/111/',dataset,'_cds_for_subset.rds'))
  target_types <-  c("CSC",
                     "LumSEC-Lip-CG",
                     "Lum-Tm4sf4",
                     "Lum-Stat4",
                     "Basal")
  
  
  
  cds_subset <- cds[,cds@colData@listData[["newcelltype"]] %in% target_types]
  #cds_subset <- choose_cells(cds)  ####手动操作####
  cds_subset <- cluster_cells(cds_subset,random_seed = 2024)
  cds_subset <- learn_graph(cds_subset,use_partition = TRUE)
  cds_subset <- order_cells(cds_subset)   ####手动操作####
  
  save_monocle_objects(cds=cds_subset, directory_path=paste0('D:/111/',dataset,'_subset_cds.rds'))
}
#######################################

#########绘图##################
if(file.exists(paste0(dataset,'_subset_cds.rds'))){
  cds_subset <- load_monocle_objects(directory_path=paste0(dataset,'_subset_cds.rds'))
  p1 <- plot_cells(cds_subset,
                   color_cells_by = "pseudotime",
                   label_cell_groups=TRUE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   graph_label_size=3)
  print(p1)
  p2 <- plot_cells(cds_subset,
                   color_cells_by = "newcelltype",
                   label_groups_by_cluster=TRUE,
                   label_leaves=TRUE,
                   label_branch_points=TRUE,
                   graph_label_size=3)
  print(p2)
  
  # 使用plot_grid来组合图表，设置一行三个图
  p1p2 <- wrap_plots(p1,p2)
  
  # 保存featureplot，调整width和height以适应新的布局
  ggsave(paste0(dataset, "_m3_pseudotime.pdf"), p1p2, width = 15, height = 10)
  
  ciliated_cds_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=cores)
  pr_deg_ids <- subset(ciliated_cds_pr_test_res, q_value < 0.05)
  
  Track_genes_sig <- pr_deg_ids %>% top_n(n=5, morans_I)%>% pull(gene_short_name) %>% as.character()
  
  
  p3 <- plot_genes_in_pseudotime(cds_subset[Track_genes_sig,],
                                 color_cells_by="newcelltype",
                                 min_expr=0.5)
  ggsave(paste0(dataset, "_m3_Jitterplot.png"), plot = p3, width = 8, height = 6)
  
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
  
  ggsave(paste0(dataset, "_m3_dotplot1.png"), plot = p4, width = 8, height = 6)
}

#####电脑上画图######
if(FALSE){
  cds <- load_monocle_objects(directory_path=paste0('D:/111/',dataset,'_subset_cds.rds'))
  target_types <-  c("MaSC",
                     "MaSC-Pro",
                     "LumSEC-Lac",
                     "LumSEC-Lip",
                     "Lum-Kit")
  cds_subset <- cds[,cds@colData@listData[["newcelltype"]] %in% target_types]
  p1 <- plot_cells(cds_subset,
                   color_cells_by = "pseudotime",
                   label_cell_groups=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=FALSE,
                   graph_label_size=3)
  print(p1)
  p2 <- plot_cells(cds_subset,
                   color_cells_by = "pseudotime",show_trajectory_graph = FALSE)
  print(p2)
  
  # 使用plot_grid来组合图表，设置一行三个图
  p1p2 <- wrap_plots(p1,p2)
  
  # 保存featureplot，调整width和height以适应新的布局
  ggsave(paste0(dataset, "_m3_pseudotime.pdf"), p1p2, width = 10, height = 5)
  ggsave(paste0(dataset, "_m3_pseudotime.png"), p1p2, width = 10, height = 5)
  
  #####################
  # ======== Step 1: 筛选你感兴趣的基因 ========
  
  genes_of_interest <- c("Dsg3")
  
  # 从 cds 中提取这些基因的 subset（注意必须是 gene_short_name 匹配）
  cds_subset <- cds_subset[rowData(cds_subset)$gene_short_name %in% genes_of_interest, ]
  
  # ======== Step 2: 绘图 ========
  p <- plot_genes_in_pseudotime(
    cds_subset,
    color_cells_by = "pseudotime",    # 可改为 "cluster" / "cell_type" / "state"
    min_expr = 0.1,                     # 可选，默认全部显示
    ncol = 1,                         # 每页显示几列子图（改为 3 横排显示）
    label_by_short_name = TRUE,
  )
  
  # ======== Step 3: 保存为 PDF 文件 ========
  ggsave(paste0(dataset,"_pseudotime_expression.pdf"), 
         plot = p,
         width = 3, height = 3, dpi = 300)  # 高度根据基因数调整
}


