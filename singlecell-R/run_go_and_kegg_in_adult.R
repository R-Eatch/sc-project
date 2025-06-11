library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(homologene)
library(patchwork)
library(tidyr)
#######################

df <- read.csv('D:/111/DEG_Analysis_Output/Stage_Intersections\\intersection_stage_3-ges-la_all.csv')
genelist <- df[df$Combination == "M-MG&R-MG&S-MG", c("Gene")]

##############################
ont <- "BP"
output_go <- TRUE
output_kegg <- TRUE
dataset <- 'ALL_MG'
group_name <- 'ges-la'
#############################


##############################
folder <- "D:/111/DEG_Analysis_Output/GO-KEGG-output/"

# 检查文件夹是否存在
if (!dir.exists(folder)) {
  # 不存在就创建，recursive = TRUE 会递归创建多级目录
  dir.create(folder, recursive = TRUE)
}
p1_list <- list()
# GO富集分析
if (output_go) {
  ea <- enrichGO(genelist,
                 OrgDb = org.Mm.eg.db,
                 keyType = "SYMBOL",
                 ont = ont,
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2)
  
  p1 <- dotplot(ea, showCategory = 15) + 
    ggtitle(paste0(dataset, '-', group_name, "_GO_Dotplot"))
  write.csv(ea@result,paste0(folder,'GO-result-',dataset,'-',group_name,'.csv'),row.names = FALSE)
  print(paste0(dataset, '-', group_name, "_GO_Dotplot finished"))
  
  # 将生成的GO图添加到p1_list中
  p1_list[[group_name]] <- p1
}
p2_list <- list()
# KEGG富集分析
if (output_kegg) {
  # 将基因符号转换为ENTREZ ID
  entrez_ids <- bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  
  kegg_result <- enrichKEGG(entrez_ids$ENTREZID,
                            organism = 'mmu',  # 小鼠
                            keyType = "ncbi-geneid",
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.1)
  p2 <- dotplot(kegg_result, showCategory = 15) + 
    ggtitle(paste0(dataset, '-', group_name, "_KEGG_Dotplot"))
  print(paste0(dataset, '-', group_name, "_KEGG_Dotplot finished"))
  entrez_ids$ENTREZID <- as.character(entrez_ids$ENTREZID)
  
  # 替换 geneID
  kegg_result@result$geneID <- sapply(kegg_result@result$geneID, function(row) {
    gene_ids <- strsplit(row, "/")[[1]]
    symbols <- sapply(gene_ids, function(id) {
      # 通过 match 查找索引
      idx <- match(id, entrez_ids$ENTREZID)
      if (!is.na(idx)) entrez_ids$SYMBOL[idx] else id
    })
    paste(symbols, collapse = "/")
  })
  write.csv(kegg_result@result,paste0(folder,'KEGG-result-',dataset,'-',group_name,'.csv'),row.names = FALSE)
  # 将生成的KEGG图添加到p2_list中
  p2_list[[group_name]] <- p2
  
}

if (output_go){
  p3=wrap_plots(p1_list)
  ggsave(paste0(dataset,"_GO_Dotplot.png"),path = folder,plot=p3,height = 30,width = 30,limitsize = FALSE)
  ggsave(paste0(dataset,"_GO_Dotplot.pdf"),path = folder,plot=p3,height = 30,width = 30,limitsize = FALSE)
}
if (output_kegg){
  p4=wrap_plots(p2_list)
  ggsave(paste0(dataset,"_KEGG_Dotplot.png"),path = folder,plot=p4,height = 30,width = 30,limitsize = FALSE)
  ggsave(paste0(dataset,"_KEGG_Dotplot.pdf"),path = folder,plot=p3,height = 30,width = 30,limitsize = FALSE)
}

