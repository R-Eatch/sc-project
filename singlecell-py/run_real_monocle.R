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
dataset <- "M-MG"
cores <- 8
s.genes <- c("Prlr",'Esr1',"Epcam",'Acta2','Lalba','Elf5','Esr1','Lef1','Ppl','Mif')
brunch_gene <- s.genes
subset_celltype <- c("ASC",'LumSEC-Lip',"LumSEC-MG-like")
root_state <- 1 #for ordercells
subset_clusters <- c()
branch_point <- 1 ##for BEAM analysis
set.seed(2024)
subset_cells_number<- 5000
pse_gene_num <- 1000
n_hvgs <- 5000
use_ordergene <- "~stage"
use_hvg=TRUE
use_seurat_hvgs = TRUE
ifsubsetcelltype <- FALSE
ifselect_number <- FALSE
color_dict <- c(
  "ASC" = "#e31a1c",           # Strong Red
  "ASC-t2-rb" = "#a6cee3",     # Light Blue
  "ASC-t3-rb" = "#fb9a99",     # Light Red/Pink
  "Basal" = "#1f78b4",         # Strong Blue
  "CSC" = "#fdbf6f",           # Light Orange/Peach
  "Epi-Fibro" = "#cab2d6",     # Lavender
  "Epi-Pro" = "#D3D3D3",       # LightGray
  "Epi-Sbsn" = "#808080",       # Gray
  "Lum-Basal" = "#778899",     # LightSlateGray
  "Lum-Immune" = "#ffd700",     # Gold
  "Lum-Kit" = "#008080",       # Teal
  "Lum-Krt8" = "#b2df8a",     # Light Green
  "Lum-Ppl" = "#000000",       # Black
  "Lum-Stat4" = "#8FBC8F",     # DarkSeaGreen (New)
  "Lum-Tm4sf4" = "#00FFFF",     # Cyan/Aqua (New)
  "LumHR" = "#ff7f00",         # Strong Orange
  "LumSEC-Hmgcr" = "#2F4F4F",  # DarkSlateGray
  "LumSEC-Lac" = "#6a3d9a",     # Strong Purple
  "LumSEC-Lip" = "#b15928",     # Brown
  "LumSEC-Lip-CG" = "#D2B48C", # Tan (New)
  "LumSEC-Lip-E" = "#4B0082",  # Indigo
  "LumSEC-Lip-I" = "#87CEEB",  # SkyBlue
  "LumSEC-MG-like" = "#B8860B",# DarkGoldenrod
  "LumSEC-Mgp" = "#00FF7F",     # SpringGreen
  "LumSEC-Vcam1" = "#FF00FF",  # Magenta/Fuchsia (Changed)
  "MaSC" = "#33a02c",         # Strong Green
  "MaSC-Pro" = "#FF69B4",     # HotPink
  "MaSC-t2-sg" = "#FAEBD7"      # AntiqueWhite
)
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

# 新函数：随机抽取指定数量的细胞并生成新Seurat对象
random_subset_cells <- function(seurat_obj, total_cells) {
  # 检查 Seurat 对象是否为空
  if (is.null(seurat_obj)) {
    stop("The Seurat object is NULL. Please provide a valid Seurat object.")
  }
  
  # 获取 Seurat 对象中所有细胞的名称
  all_cells <- Cells(seurat_obj)
  
  # 检查 Seurat 对象中的细胞数量是否足够
  if (length(all_cells) < total_cells) {
    stop("The total number of cells in the Seurat object is less than the specified number of cells to be sampled.")
  }
  
  # 随机抽取指定数量的细胞
  selected_cells <- sample(all_cells, total_cells)
  
  # 根据选出的细胞创建一个新的 Seurat 对象
  final_subset <- subset(seurat_obj, cells = selected_cells)
  
  # 返回最终的 Seurat 对象
  return(final_subset)
}
# 函数定义
filter_color_scheme <- function(color_dict, seurat_object) {
  existing_celltypes <- unique(seurat_object$newcelltype)

  filtered_colors <- color_dict[names(color_dict) %in% existing_celltypes]
  
  # 返回过滤后的配色向量
  return(filtered_colors)
}

##################


ea <- readRDS(paste0(dataset,"_for_monocle.rds"))
ea$celltype=ea$newcelltype
if(ifsubsetcelltype){
  ea <- subset(ea, subset = celltype %in% subset_celltype)
}else{
  subset_celltype <- unique(ea$newcelltype)
}
#ea <- random_subset_cells(ea,subset_cells_number)
filtered_colors <- filter_color_scheme(color_dict, ea)
if(ifselect_number){
  ea <- random_subset_by_celltype(ea,subset_celltype,subset_cells_number)
}
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
if (use_hvg){
  if(use_seurat_hvgs){
    ea <- NormalizeData(ea, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    ea <- FindVariableFeatures(ea, selection.method = "vst", nfeatures = n_hvgs, verbose = FALSE)
    express_genes <- VariableFeatures(ea)
  }else{
    express_genes <-rownames(ea[["RNA"]]@meta.features)[ea[["RNA"]]@meta.features$highly.variable == TRUE]
  }
}else{
  express_genes <-rownames(ea[["RNA"]]@meta.features)
}
diff_test_res <- differentialGeneTest(cds[express_genes,],
                                      fullModelFormulaStr = use_ordergene,cores = cores)
ordering_genes <- row.names(subset(diff_test_res, qval < 0.001))
cds <- setOrderingFilter(cds, ordering_genes)
p11 <- plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')
cds <- orderCells(cds)
cds <- orderCells(cds,root_state = root_state)####定义根的起点
p1 <- plot_cell_trajectory(cds, color_by = "newcelltype", size = 10)+ scale_color_manual(values = filtered_colors)
p2 <- plot_cell_trajectory(cds,color_by="Pseudotime", size=10,show_backbone=TRUE) + scale_color_viridis_c(option = 'plasma')
p3 <- plot_cell_trajectory(cds, color_by = "stage") + scale_color_viridis_d(option = 'D')
p3.5 <- plot_cell_trajectory(cds, color_by = "State")
p3.6 <- plot_cell_trajectory(cds, color_by = "leiden",cell_name_size=10)
p3.7 <- plot_cell_trajectory(cds, color_by = "sample",cell_name_size=10)
p4 <- plot_genes_jitter(cds[s.genes,], grouping = "newcelltype", color_by = "newcelltype")+
  ggtitle(paste0(dataset, " Gene Jitterplot")) +
  theme(
    legend.title = element_text(size = 10),  # 调整图例标题的大小
    legend.text = element_text(size = 10),   # 调整图例文本的大小
    plot.title = element_text(size = 10, hjust = 0.5)  # 调整标题的大小并居中
  ) +
  guides(colour = guide_legend(override.aes = list(size = 5))) + scale_color_manual(values = filtered_colors)# 调整图例中点的大小
p5 <- plot_genes_violin(cds[s.genes,], grouping = "newcelltype",color_by = 'newcelltype') + scale_color_manual(values = filtered_colors)
p6 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "newcelltype") + scale_color_manual(values = filtered_colors)
plotc1 <- (p1 | p2 | p3) /
          (p3.5 | p3.6 | p3.7)
ggsave(paste0(dataset,"pseu_fig1.pdf"), plot = p1, width = 8, height = 8,dpi=300)
ggsave(paste0(dataset,"pseu_fig2.pdf"), plot = p2, width = 8, height = 8,dpi=300)
ggsave(paste0(dataset,"pseu_fig3.pdf"), plot = p3, width = 8, height = 8,dpi=300)
ggsave(paste0(dataset,"pseu_fig4.pdf"), plot = p3.5, width = 8, height = 8,dpi=300)
ggsave(paste0(dataset,"pseu_fig5.pdf"), plot = p3.6, width = 8, height = 8,dpi=300)
ggsave(paste0(dataset,"pseu_fig6.pdf"), plot = p3.7, width = 8, height = 8,dpi=300)
plotc1 <- plotc1 + 
  plot_annotation(title = paste0(dataset," Pseudotime trajectory"),
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
plotc2 <- p4|p5|p6
ggsave(paste0(dataset,"_pseudotime.png"), plot = plotc1, width = 30, height = 20)
ggsave(paste0(dataset,"_pseudotime.pdf"), plot = plotc1, width = 30, height = 20,dpi = 300)
ggsave(paste0(dataset,"_Genes_Jitterplot.png"), plot = plotc2, width = 16, height = 8)
ggsave(paste0(dataset,"_Genes_Jitterplot.pdf"), plot = plotc2, width = 16, height = 8,dpi = 300)

if(run <- TRUE){
  diff_test_res1 <- differentialGeneTest(cds[ordering_genes,],
                                         fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = cores)
  write.csv(diff_test_res1, paste0(dataset,"all_ordergene.csv"), row.names = F)                                      
  top100 <- top_n(diff_test_res1, n=pse_gene_num, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  p <- plot_pseudotime_heatmap(cds[top100,],
                               num_clusters = 5,
                               cores = cores,return_heatmap=T,
                               show_rownames = F)
  
  clusters <- cutree(p$tree_row, k = 5) 
  clustering <- data.frame(clusters) 
  clustering[,1] <- as.character(clustering[,1]) 
  colnames(clustering) <- "Gene_Clusters" 
  table(clustering)
  write.csv(clustering, paste0(dataset,"_pseudotime_cluster.csv"), row.names = T)
  top100 <- top_n(diff_test_res1, n=pse_gene_num, desc(qval)) %>% pull(gene_short_name) %>% as.character()
  p7 <- plot_pseudotime_heatmap(cds[top100,],
                                num_clusters = 5,
                                cores = cores,return_heatmap=T,
                                show_rownames = F)
  ggsave(paste0(dataset,"_Gene_heatmap.png"), p7, width = 6.5, height = 10)
  ggsave(paste0(dataset,"_Gene_heatmap.pdf"), p7, width = 6.5, height = 10,dpi = 300)
  
}
saveRDS(cds,paste0(dataset,"_monocle.rds"))
if(TRUE){
##BEAM analysis
BEAM_res <- BEAM(cds[ordering_genes,], branch_point = branch_point, cores = cores,progenitor_method = "duplicate")
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_GENE <- subset(BEAM_res, qval<0.01) %>% pull(gene_short_name) %>% as.character()

write.csv(BEAM_res, paste0(dataset,"_point",branch_point,"_BEAM_result.csv"), row.names = F)

top100 <- top_n(BEAM_res, n=100, desc(qval)) %>% pull(gene_short_name) %>% as.character()

p8 <- plot_genes_branched_heatmap(cds[top100,],
                            branch_point = branch_point,
                            num_clusters = 3,
                            cores = cores,return_heatmap = T,
                            use_gene_short_name = T,
                            show_rownames = T)
ggsave(paste0(dataset,"_point",branch_point,"_Branch_heatmap.png"), p8$ph_res, width = 6.5, height = 10)
ggsave(paste0(dataset,"_point",branch_point,"_Branch_heatmap.png"), p8$ph_res, width = 6.5, height = 10,dpi = 300)
df1 <- p8$annotation_row
df3 <- data.frame("gene" = rownames(df1), "cluster" = df1$Cluster)
write.csv(df3,paste0(dataset,"_point",branch_point,"_Branch_cluster.csv"),row.names = F)
intrest_gene <- row.names(subset(fData(cds),
                               gene_short_name %in% brunch_gene))
p9 <- plot_genes_branched_pseudotime(cds[intrest_gene,],
                               branch_point = branch_point,
                               color_by = "newcelltype",
                               ncol = 2)+
  ggtitle(paste0(dataset, " Branch Jitterplot"))+ scale_color_manual(values = filtered_colors)
ggsave(paste0(dataset,"_point",branch_point,"_Branch_Jitterplot.png"), p9, width = 6.5, height = 10)
ggsave(paste0(dataset,"_point",branch_point,"_Branch_Jitterplot.pdf"), p9, width = 6.5, height = 10,dpi = 300)
}



