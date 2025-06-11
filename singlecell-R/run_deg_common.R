library(Seurat)
library(dplyr)
library(ggplot2)
dataset <- "R-MAG"
ea <- readRDS(paste0("../",dataset,"CCA_combined.rds"))
ea <- subset(ea,celltype == "Luminal")
ea$newstage <- ifelse(ea$time %in% c(1,2,3,4, 5, 6), 'PUB',
                   ifelse(ea$time %in% c(7,8,9), 'GES-LA','path0'))
ea$newstage_gland <- paste(ea$newstage, ea$gland, sep="-")
Idents(ea) <- ea$newstage
A <- FindAllMarkers(ea,only.pos = TRUE,logfc.threshold = 0.25,min.pct = 0.25)

Idents(ea) <- ea$newstage_gland
B <- FindConservedMarkers(ea, ident.1 = "GES-LA-MG",grouping.var = 'newstage',
                                       logfc.threshold = 0.25, min.pct = 0.25,only.pos = TRUE)
B$gene <- rownames(B)
df1 <- FindAllMarkers(subset(ea,gland == "AG"),only.pos = TRUE,logfc.threshold = 0.25,min.pct = 0.25)
df2 <- FindAllMarkers(subset(ea,gland == "MG"),only.pos = TRUE,logfc.threshold = 0.25,min.pct = 0.25)
df1 <- filter(df1,cluster == 'GES-LA-AG')
df2 <- filter(df2,cluster == 'GES-LA-MG')
A <- filter(A,p_val_adj < 0.01)
B <- filter(B,`GES-LA_p_val_adj` < 0.01)
df3 <- merge(df1,df2,by = "gene")
C <-merge(A,B,by = "gene")
E <- merge(df3,C,by = "gene")

top10_genes <- E %>%
  arrange(p_val_adj) %>%
  slice(1:20)
genelist <- top10_genes$gene

p1 <- DotPlot(ea,features = top10_genes$gene,scale = TRUE,group.by = "newstage")
p2 <- DotPlot(ea,features = top10_genes$gene,scale = TRUE,group.by = "newstage_gland")

ggsave(paste0(dataset,"newstage_DotPlot.png"), plot = p1, width = length(genelist1), height = 8, dpi = 300,bg = 'white')
ggsave(paste0(dataset,"newstage_DotPlot.png"), plot = p2, width = length(genelist1), height = 8, dpi = 300,bg = 'white')
