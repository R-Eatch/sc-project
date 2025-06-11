if (TRUE){
library(Seurat)
library(dplyr)
library(ggplot2)
library(wesanderson)
library(viridis)
dataset <- 'MRS-MG'
deg_list <- list()
FP <- wes_palette("Zissou1", 5, type = "discrete")
ea <- readRDS(paste0("../",dataset,"CCA_combined.rds"))
ea <- subset(ea,subset = time %in% c(7,8,9))
Idents(ea) <- ea$celltype
for (species in unique(ea$species)){
  seurat_obj <- subset(ea ,subset = species == species)
  df <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  df1 <- df %>%
    filter(cluster == 'Luminal' & p_val < 0.05)
  deg_list[[species]] <- df1
}
conserved_genes <- Reduce(intersect, lapply(deg_list, function(x) x$gene))
ea <- subset(ea,idents = 'Luminal')
Idents(ea) <- ea$species
DefaultAssay(ea) <- "RNA"

markers <- FindAllMarkers(ea, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
unique <- unique(as.character(top5_markers$gene))
  p_dot1 <- DotPlot(ea, features = unique) + RotatedAxis() +
    scale_colour_gradientn(colors = FP) +
    ggtitle('unique genes between species in MG' )
  p_dot2 <- DotPlot(ea, features = conserved_genes) + RotatedAxis() +
    scale_colour_gradientn(colors = FP) +
    ggtitle('common genes between species in MG')
ggsave(paste0(dataset, "species_unique_DotPlot.png"), plot = p_dot1, width = 15, height = 8)
ggsave(paste0(dataset, "species_common_DotPlot.png"), plot = p_dot2, width = 15, height = 8)
write.csv(markers, paste0(dataset, "unique_DEG.csv"))
write.csv(conserved_genes, paste0(dataset, "common_DEG.csv"))

p1 <- DimPlot(ea, reduction = "umap", group.by = "merge_dataset") +
  ggtitle(paste0(dataset, " - Subset - By Dataset"))

p2 <- DimPlot(ea, reduction = "umap", group.by = "celltype", label = TRUE, pt.size = 1) + ggtitle(paste0(dataset," - Subset - By CellType"))

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
ggsave(paste0(dataset,"_Subset_celltype.png"), plot = plotc2, width = 30, height = 20)
}



library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
######################
dataset <- 'MRS-MG'
use_select_gene <- TRUE
genelist <- c("Bmp4", "Pthlh", "Dhh", "Gli3","Clu","Slc28a3",'Irak2','Plod2','Rarb','Fabp5','Tnfrsf19','Prrx1','Cntn3','Hspg2','Ctnnal1','Lgr5','Cap2','Fam20a')
##################################
plot_all1 <- function(genelist,ea_subset,dataset){
  DefaultAssay(ea_subset) <- "RNA"
  p1 <- DotPlot(ea_subset, features = genelist, group.by = "stage_species") +
    ggtitle(paste0("Expression patterns across stages and glands in ",dataset)) +
    theme_minimal() +
    scale_color_gradient(low = "green", high = "red") +  
    scale_size(range = c(3, 10))
  all_plots <- lapply(genelist, function(gene) {
    data_to_plot <- FetchData(ea_subset, vars = c("stage", "species", gene),layer = "data")
    
    # 处理数据：按时间和腺体类型分组计算平均表达
    average_expression <- data_to_plot %>%
      group_by(stage, species) %>%
      summarise(Avg_Expression = mean(!!sym(gene), na.rm = TRUE), .groups = 'drop')
    plot <- ggplot(average_expression, aes(x = stage, y = Avg_Expression, color = species, group = species)) +
      geom_line() +
      geom_point() +
      labs(title = paste(gene),
           subtitle = "Expression in different species types",
           x = "Time",
           y = "Average Expression") +
      theme_minimal() +
      theme(text = element_text(size = 10),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "top")
    
    return(plot)
  })
  p2 <- wrap_plots(all_plots, ncol = 4,nrow = (length(genelist)%/%4)+1)
  
  all_plots <- lapply(genelist, function(gene) {
    data_to_plot <- FetchData(ea_subset, vars = c("time", "species", gene),layer = "data")
    
    # 处理数据：按时间和腺体类型分组计算平均表达
    average_expression <- data_to_plot %>%
      group_by(time, species) %>%
      summarise(Avg_Expression = mean(!!sym(gene), na.rm = TRUE), .groups = 'drop')
    # 绘制折线图，为每种腺体类型绘制一条线
    plot <- ggplot(average_expression, aes(x = time, y = Avg_Expression, color = species, group = species)) +
      geom_line() +
      geom_point() +
      labs(title = paste(gene),
           subtitle = "Expression in different species types",
           x = "Time",
           y = "Average Expression") +
      theme_minimal() +
      theme(text = element_text(size = 10),
            plot.title = element_text(face = "bold", hjust = 0.5),
            legend.position = "top")
    
    return(plot)
  })
  p3 <- wrap_plots(all_plots, ncol = 4)
  return(list(p1=p1,p2=p2,p3=p3))
}
##############################
ea <- readRDS(paste0("../",dataset,"CCA_combined.rds"))
Idents(ea) <- ea$celltype
#ea <- subset(ea,idents = "Luminal")
stage1.name <- 'Stage1'
stage2.name <- 'Stage2'
stage3.name <- 'Stage3'
ea$stage <- ifelse(ea$time %in% c(4, 5,6), 'Stage1',
                   ifelse(ea$time %in% c(7, 8), 'Stage2', ifelse(ea$time %in% c(9), 'Stage3', 'Stage0')))


ea$stage_species <- paste(ea$stage, ea$species, sep="-")

Idents(ea) <- ea$stage_species

############### select-gene ######################################
if (use_select_gene){
  result <- plot_all1(genelist,ea,dataset)
  ggsave(paste0(dataset,"select_Genes_DotPlot.png"), plot = result$p1, width = length(genelist), height = 8, dpi = 300,bg = 'white')
  ggsave(paste0(dataset,"select_expression_plots_stage_series.png"), plot = result$p2, width = 16, height = ((length(genelist)%/%4)+1)*4)
  ggsave(paste0(dataset,"select_expression_plots_time_series.png"), plot = result$p3, width = 16, height = ((length(genelist)%/%4)+1)*4)
}
####################DEG-gene##################################
genelist1 <- markers %>%
  group_by(cluster) %>%
  top_n(n=20)
  result <- plot_all1(genelist1$gene,ea,dataset)
  ggsave(paste0(dataset,"DEG_Genes_DotPlot.png"), plot = result$p1, width = length(genelist), height = 8, dpi = 300,bg = 'white')
  ggsave(paste0(dataset,"DEG_expression_plots_stage_series.png"), plot = result$p2, width = 16, height = ((length(genelist)%/%4)+1)*4)
  ggsave(paste0(dataset,"DEG_expression_plots_time_series.png"), plot = result$p3, width = 16, height = ((length(genelist)%/%4)+1)*4)

