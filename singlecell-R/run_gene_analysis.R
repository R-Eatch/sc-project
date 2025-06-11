library(dplyr)
library(readr)
library(homologene)
library(purrr)
library(tidyr)
###########################
datasets <- c("M_MG","R_MG","S_MG")
filename <- "BL_MG_MRS"
all_samples <- tibble()
############################
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
#############################

for (name in datasets){
  data <- read_csv(paste0("C:/Users/Lenovo/Desktop/",name,"Branch_all.csv"))
  if (name == "M_MG"){
    data <- convert_mouse_to_human_expression(data)
  }
  data$dataset <- name
  all_samples <- bind_rows(all_samples, data)
}
gene_stats <- all_samples %>%
  group_by(gene_short_name) %>%
  summarise(
    combined_pval = -2 * sum(log(pval), na.rm = TRUE),
    datasets = list(unique(dataset)),  # 保留所有包含该gene的dataset
    count = n_distinct(dataset),       # 计算每个基因出现在不同样本的次数
    .groups = 'drop'
  ) %>%
  mutate(
    fisher_pvalue = pchisq(combined_pval, df = 2 * count, lower.tail = FALSE)  # 使用正确的自由度计算 Fisher p-value
  )

# 筛选出所有样本中都有的基因
num <- length(datasets)
common_genes <- gene_stats %>%
  filter(count == num) %>%
  select(gene_short_name, fisher_pvalue,datasets) %>%
  arrange(fisher_pvalue)


# 筛选出只在一个样本中出现的基因，保留该gene所在的dataset
unique_genes <- gene_stats %>%
  filter(count == 1) %>%
  mutate(datasets = sapply(datasets, toString)) %>%  # 确保 datasets 是字符型
  select(genename = gene_short_name, fisher_pvalue, datasets)

write_csv(common_genes, paste0(filename,"common_genes.csv"))
write_csv(unique_genes, paste0(filename,"unique_genes.csv"))
