library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(homologene)
library(patchwork)
library(tidyr)
#######################
convert_mouse_to_human_expression <- function(mouse_expression) {
  mouse_genes <- mouse_expression$gene_short_name
  new_gene_names <- mouse_genes
  # 使用 mouse2human 函数进行基因名转换
  human_genes <- mouse2human(genes = mouse_genes, db = homologeneData)
  human_genes_unique_first_pass <- human_genes[!duplicated(human_genes$humanGene), ]
  human_genes_unique <- human_genes_unique_first_pass[!duplicated(human_genes_unique_first_pass$mouseGene), ]
  new_gene_names <- rep(NA, length(mouse_genes))
  
  # 遍历所有小鼠基因，查找并尝试替换为对应的人类基因名
  for (i in seq_along(mouse_genes)) {
    if (mouse_genes[i] %in% human_genes_unique$mouseGene) {
      new_gene_names[i] <- human_genes_unique$humanGene[human_genes_unique$mouseGene == mouse_genes[i]]
    } else {
      # 如果找不到对应的人类基因，保留原始的小鼠基因名
      new_gene_names[i] <- mouse_genes[i]
    }
  }
  mouse_expression$gene_short_name <- new_gene_names
  
  return(mouse_expression)
}
#########################
datasetlist <-list('M-MG','R-MG','R-AG','R-CG','S-MG','S-AG')
for (dataset in datasetlist){
  data <- read.csv(paste0("D:/111/DEG-result/",dataset,"_ranked_genes.csv"))
  ont <- "BP"
  isMOUSE <- FALSE
  output_go <- TRUE
  output_kegg <- FALSE
  n_gene <- 999999
  colnames(data)[4] <- "pvalue"  # 确保第二列是pvalue列
  groups <- unique(data$group)  # 获取所有唯一的组别
  genelists <- list()  # 创建一个空列表存储结果
  p1_list <- list()  # 用于GO富集图
  p2_list <- list()  # 用于KEGG富集图
  for (group_name in groups) {
    # 对每个组的数据进行筛选
    group_data <- data %>% filter(group == group_name)
    # 进一步筛选pvalue < 0.01的基因
    group_data <- group_data %>% filter(pvalue < 0.01)
    # 排序pvalue，并选择前100行
    go1 <- group_data %>% arrange(desc(logfoldchange)) %>% head(n_gene)  # 按pvalue升序排列，取前100行
    
    # 提取基因名，并存储到genelist中
    genelist <- go1[,c("gene","logfoldchange")]
    gene_scores <- genelist[[2]]             
    names(gene_scores) <- genelist[[1]] 
    # GO富集分析
    if (output_go) {
      ea <- gseGO(gene_scores,
                     OrgDb = org.Mm.eg.db,
                     keyType = "SYMBOL",
                     ont = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05)
      
      p1 <- dotplot(ea, showCategory = 15) + 
        ggtitle(paste0(dataset, '-', group_name, "_GO_Dotplot"))
      write.csv(ea@result,paste0('D:/111/DEG-result/','GO-result-',dataset,'-',group_name,'.csv'),row.names = FALSE)
      print(paste0(dataset, '-', group_name, "_GO_Dotplot finished"))

      # 将生成的GO图添加到p1_list中
      p1_list[[group_name]] <- p1
    }
    
    # KEGG富集分析
    if (output_kegg) {
      # 将基因符号转换为ENTREZ ID
      entrez_ids <- bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
      
      kegg_result <- gseKEGG(entrez_ids$ENTREZID,
                                organism = 'mmu',  # 小鼠
                                keyType = "ncbi-geneid",
                                pAdjustMethod = "BH",
                                pvalueCutoff = 0.05)
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
      write.csv(kegg_result@result,paste0('D:/111/DEG-result/','KEGG-result-',dataset,'-',group_name,'.csv'),row.names = FALSE)
      # 将生成的KEGG图添加到p2_list中
      p2_list[[group_name]] <- p2
      
    }
  }
  if (output_go){
    p3=wrap_plots(p1_list)
    ggsave(paste0(dataset,"_GO_Dotplot.png"),path = "D:/111/DEG-result",plot=p3,height = 30,width = 30,limitsize = FALSE)
  }
  if (output_kegg){
    p4=wrap_plots(p2_list)
    ggsave(paste0(dataset,"_KEGG_Dotplot.png"),path = "D:/111/DEG-result",plot=p4,height = 30,width = 30,limitsize = FALSE)
  }
}
