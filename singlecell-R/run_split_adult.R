# 加载Seurat包
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
datasetlist <- c('M_MG','R_MG','S_MG','R_AG','R_CG','S_AG')
cell_percent_all <- data.frame()
plotlist <- list()

lapply(datasetlist,function(dataset) {
  if (dataset %in% c('S_MG','S_AG')){
    dataset1 <- gsub("_","-",dataset)
    path <- paste0('../',dataset1,'/1.subset/',dataset,'subset_for_monocle_CCA.rds')
  }else{
    path <- paste0('../',dataset1,'/1.subset/',dataset,'subset_for_monocle.rds') 
  }
  print(paste0('cruuent dataset : ',dataset))
  ea <- readRDS(path)
  ea_subset <- subset(ea,subset = time == '6')
  p1 <- DimPlot(ea_subset,group.by = 'celltype') +
    ggtitle(paste0(dataset,'UMAP plot in adult(time6)'))
  Idents(ea_subset) <- ea_subset$celltype
  markers <- FindAllMarkers(ea_subset,logfc.threshold = 0.25,min.pct = 0.25,only.pos = TRUE)
  top5 <- markers %>% 
    group_by(cluster) %>%
    slice_head(n=5)
  
  p2 <- DotPlot(ea_subset,group.by = 'celltype',features = top5$gene) + 
    ggtitle(paste0(dataset,' Dotplot in adult(time6)')) +
    scale_color_viridis(option = "turbo") + 
    scale_fill_viridis(option = "turbo") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  cell_number <- data.frame(table(ea_subset$celltype))
  cell_percent <- data.frame(prop.table((table(ea_subset$celltype))))
  colnames(cell_percent) <- c('celltype','percent')
  colnames(cell_number) <- c('celltype','number')
  df1 <- merge(cell_percent,cell_number)
  df1$name <- dataset
  print(path)
  plotlist[[dataset]] <- list(dimplot = p1,dotplot = p2)
  cell_percent_all <- rbind(cell_percent_all, df1)
  }
)
for (dataset in datasetlist){
  ggsave(plot = plotlist[[dataset]]$dotplot,filename = paste0(dataset,'Dotplot.png'),width = 30,height = 20)
}
dim_plots <- lapply(plotlist, function(x) x$dimplot)
patchwork::wrap_plots(dim_plots, ncol = 2) + plot_annotation(tag_levels = 'A')
ggsave("umap_plots.png", width = 15, height = 10)

bar_plot <-   ggplot(cell_percent_all, aes(x = name, y = percent, fill = celltype)) +
      geom_col(position = "fill") + 
      facet_wrap(~ name, scales = "free_x") + 
      labs(x = "dataset", y = "percent",title = "Cell Proportion in Adult Samples") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = bar_plot,filename = 'barplot.png',width = 30,height = 20)
write.csv(cell_percent_all,file = 'cell_percent.csv')


