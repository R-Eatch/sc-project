library(Seurat)
library(ggplot2)
library(wesanderson)
library(patchwork)
library(homologene)
###########################
featureall <- function(seurat_object, genelist, dataset_name) {
  plots <- list()  # 创建一个空列表来存储图形对象
  
  # 首先添加一个基于 CellType 的 UMAP 绘图，标题为数据集名
  if ("umap" %in% names(seurat_object@reductions)) {
    umap_plot <- DimPlot(seurat_object, group.by = "celltype", label = TRUE, pt.size = 1) +
      ggtitle(paste("UMAP of", dataset_name))  # 设置基于数据集名称的标题
    plots[["UMAP"]] <- umap_plot
  } else {
    print("UMAP数据未找到，请先运行降维。")
  }
  
  # 绘制每个基因的特征图
  for (gene in genelist) {
    if (gene %in% rownames(seurat_object@assays[["RNA"]]@counts)) {
      p1 <- FeaturePlot(seurat_object, features = gene, cols = FP, pt.size = 1, min.cutoff = "q10", max.cutoff = "q90")
      print(paste0("有:", gene))
      plots[[gene]] <- p1  # 存储图形对象到列表
    } else {
      print(paste0("没有:", gene))
    }
  }
  
  # 如果有生成的图形，组合并保存它们
  if (length(plots) > 1) {  # 确保至少有两个图，包括UMAP图
    # 每行显示5个图
    n_col <- 5
    n_row <- ceiling(length(plots) / n_col)  # 计算需要的行数
    
    # 计算总体高度和宽度
    height_per_plot <- 5  # 每个图的高度为5英寸
    width_per_plot <- 7   # 每个图的宽度为7英寸
    total_height <- n_row * height_per_plot
    total_width = min(n_col, length(plots)) * width_per_plot  # 总宽度取决于列数或图形的数量
    
    # 使用 patchwork 的 wrap_plots 和 plot_layout 来组合图形
    combined_plot <- wrap_plots(plots) + 
      plot_layout(nrow = n_row, ncol = n_col)
    
    # 使用 ggsave 保存图形，确保每个子图有足够的空间
    file_name <- paste("combined_feature_plots_", dataset_name, ".png", sep = "")
    ggsave(plot = combined_plot, filename = file_name, width = total_width, height = total_height, limitsize = FALSE)
    print(paste("图形已保存为", file_name))
  } else {
    print("没有可显示的图形。")
  }
}

###########################
ea <- readRDS("D:/111/S-AGsubset_for_monocle.rds")
dataset <- "S_AG"
isMOUSE <- TRUE
getwd()
genelist <- c("CDK1","TOP2A","MKI67","SLPI","KRT18","ERBB4","KRT15","ACTA2","ESR1","STAT5","KRT1","PRLR","KRT5","EPCAM","ACTG2","ACTA1","KRT14","ACTG1","FABP4","ELOVL3")
FP <- wes_palette("Zissou1", 5, type = "discrete")
genelist <- c("CADPS2", "CLU", "ERBB4", "FABP5", "KRT5", "LMO7", "MAST4", "MGST1", 
              "PRKG1", "PTPRK", "THBS1", "ALCAM", "TBC1D9", "PGR", "KLF5", "SORBS2", 
              "GATA3", "PDGFC", "F3", "KCNMA1", "RASEF", "EPCAM", "CDH1", "CD36", 
              "KRT7", "SLC1A3", "EHF", "SEMA3C", "AR", "CERS6", "ANGPT1", "ENPP3", 
              "BARX2", "NCKAP5", "PDE7B", "NRXN1", "GCLC", "GLIS3", "TGFB2", "DSC3", 
              "IGFBP5", "WWC1", "PTPRM", "CRYAB", "CHL1", "DPYD", "RSAD2", "DCLK1", 
              "HOMER1", "CTSC", "SOX6", "LGMN", "KIF26B", "SLC7A11", "FRMD4A", "PREX2", 
              "NFATC2", "CD74", "LAMA3", "POSTN", "MYH11", "ADAMTS5")##HS_MG_AG
genelist <- c("ARHGAP20", "CACNA1C", "CADPS2", "CLU", "COL12A1", "COL3A1", 
              "CRYAB", "DGKI", "EPCAM", "ERBB4", "FABP5", "FBN1", "IGFBP5", 
              "IGFBP6", "KCNMA1", "KRT17", "KRT5", "KRT7", "KRT8", "LAMA3", 
              "LMOD1", "LPL", "MEF2C", "MGST1", "MME", "MYH11", "NFATC2", 
              "PDGFC", "PGM5", "POSTN", "PRKG1", "PTN", "RASEF", "RHOJ", 
              "RRAD", "SEMA5A", "SLC1A3", "SYNPO2", "THBS1", "TNC")##BL_MG_AG
genelist <- unique(c("FABP5", "DSC3", "CDK1", "TOP2A", 
                     "KRT1", "COL3A1", "FBN1", "CLU", 
                     "CADPS2", "SOX6", "ACTA2", "PRKG1", 
                     "MYH11", "KRT5"))
##interest in BL_MG_AG
genelist <- c("ACTG2", "CACNA1C", "CAV1", "CDH13", "CLU", "COL12A1", 
              "COL3A1", "COL4A1", "COL4A2", "DCN", "ERBB4", "FABP3", 
              "FHOD3", "HSPB1", "IGFBP5", "IGFBP6", "KCNMA1", "KRT17", 
              "KRT5", "LAMA1", "LAMA3", "LMOD1", "LYZ", "MME", "MYH11", 
              "NFATC2", "PALLD", "PCDH7", "PGM5", "PLIN2", "POSTN", 
              "PRKG1", "PTPRZ1", "RASEF", "RRAD", "SEMA5A", "SLC1A3", 
              "SLC35F1", "SYNPO2", "XDH", "CSRP1", "ANTXR1", "EPCAM", 
              "KITLG")#


if(isMOUSE){
  humangenelist <- human2mouse(genelist)
  genelist <- humangenelist$mouseGene
} else{
  genelist <- genelist
}

#############################
featureall(seurat_object = ea, genelist = genelist , dataset_name = dataset)


