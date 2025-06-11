# 清理环境（可选）
 rm(list = ls()) 


library(org.Mm.eg.db)
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse) 
library(stringr)
library(ggprism) 

# --- (gground is not needed) ---

#########################
# 2. 定义数据集列表和参数 
datasetlist <- list('M-MG',
                    'R-MG',
                    'R-AG',
                    'R-CG',
                    'S-MG',
                    'S-AG')
n_gene <- 500 
n_top_pathways <- 5 # 每个类别取 top N 通路
output_dir <- "D:/111/DEG-result/" 
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 定义颜色方案 (Requirement 4)
pal <- c('#eaa052', '#b74147', '#90ad5b', '#23929c')
# 确保名称与固定顺序 'BP', 'CC', 'MF', 'KEGG' 对应
names(pal) <- c('BP', 'CC', 'MF', 'KEGG') 

# 定义左侧标签条的固定宽度和间距单位 (Requirement 1b)
label_bar_width_unit <- 1.0 # 控制X轴方向宽度 (-3*width to -2*width)
label_bar_x_start <- -3 * label_bar_width_unit
label_bar_x_end <- -2 * label_bar_width_unit
gene_count_x_pos <- -1 * label_bar_width_unit # X位置放基因count圆圈

# 3. 开始循环处理每个数据集 
for (dataset in datasetlist) {
  print(paste("Processing dataset:", dataset))
  data_path <- file.path("D:/111/DEG-result/", paste0(dataset, "_ranked_genes.csv"))
  
  if (!file.exists(data_path)) {
    print(paste("Skipping dataset:", dataset, "- File not found:", data_path))
    next 
  }
  
  data <- read.csv(data_path)
  names(data)[names(data)=='pval'] <- 'pvalue'
  
  if (!all(c("group", "gene", "logfoldchange", "pvalue") %in% colnames(data))) {
    print(paste("Error: Required columns not found in", dataset))
    next
  }
  
  groups <- unique(data$group) 
  
  # 4. 循环处理每个细胞类型 (group) 
  for (group_name in groups) {
    print(paste("  Processing group:", group_name))
    
    # 筛选上调基因
    group_data <- data %>% filter(group == group_name)
    go1 <- group_data %>% 
      filter(logfoldchange > 0) %>% 
      arrange(pvalue) %>% 
      head(n_gene)
    
    genelist <- as.character(go1$gene)
    
    if (length(genelist) == 0) {
      print(paste("    No upregulated genes found for", group_name, "- Skipping enrichment."))
      next 
    }
    
    # 5. GO富集分析 (ont = "ALL")
    print("    Running GO enrichment...")
    ea_all <- NULL
    tryCatch({
      ea_all <- enrichGO(gene=genelist, OrgDb=org.Mm.eg.db, keyType="SYMBOL", ont="ALL", 
                         pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE) 
    }, error = function(e) { print(paste("    Error GO:", e$message)) })
    
    # 6. KEGG富集分析 
    print("    Running KEGG enrichment...")
    kegg_result <- NULL 
    entrez_ids_map <- bitr(genelist, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db, drop = TRUE)
    
    if (nrow(entrez_ids_map) > 0) {
      tryCatch({
        kegg_result_raw <- enrichKEGG(gene=entrez_ids_map$ENTREZID, organism='mmu', keyType="ncbi-geneid", 
                                      pAdjustMethod="BH", pvalueCutoff=0.1, qvalueCutoff=0.2)
        if (!is.null(kegg_result_raw) && nrow(kegg_result_raw) > 0) {
          kegg_result <- setReadable(kegg_result_raw, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
        }
      }, error = function(e) { print(paste("    Error KEGG:", e$message)) })
    } else {
      print("    No genes mapped to Entrez IDs for KEGG.")
    }
    
    # 7. 数据处理与绘图准备 
    print("    Preparing data for plot...")
    
    # 转数据框 & 保存原始结果
    GO_df <- NULL
    if (!is.null(ea_all) && nrow(ea_all@result) > 0) {
      GO_df <- as.data.frame(ea_all)
      write.csv(GO_df, file.path(output_dir, paste0('GO-result-', dataset, '-', group_name, '.csv')), row.names = FALSE)
    } else { print("    No significant GO terms.") }
    
    KEGG_df <- NULL
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      KEGG_df <- as.data.frame(kegg_result)
      # 清理 KEGG 描述
      KEGG_df$Description <- gsub("\\s*-\\s*Mus musculus.*$", "", KEGG_df$Description)
      KEGG_df$Description <- gsub("\\s*\\(house mouse\\)$", "", KEGG_df$Description)
      KEGG_df$Description <- trimws(KEGG_df$Description) 
      write.csv(KEGG_df, file.path(output_dir, paste0('KEGG-result-', dataset, '-', group_name, '_cleaned.csv')), row.names = FALSE)
    } else { print("    No significant KEGG terms.") }
    
    # 筛选 TOP 通路
    use_pathway_list <- list()
    if (!is.null(GO_df)) {
      go_top <- GO_df %>% filter(ONTOLOGY %in% c("BP", "CC", "MF")) %>% 
        group_by(ONTOLOGY) %>% arrange(p.adjust, desc(Count)) %>% 
        slice_head(n = n_top_pathways) %>% ungroup()
      if(nrow(go_top) > 0) use_pathway_list[['GO']] <- go_top
    }
    if (!is.null(KEGG_df)) { 
      kegg_top <- KEGG_df %>% arrange(p.adjust, desc(Count)) %>% 
        slice_head(n = n_top_pathways) %>% mutate(ONTOLOGY = 'KEGG') 
      if(nrow(kegg_top) > 0) use_pathway_list[['KEGG']] <- kegg_top
    }
    
    # 合并并排序 (Requirement 2: BP->CC->MF->KEGG top-to-bottom)
    if (length(use_pathway_list) > 0) {
      use_pathway_group <- bind_rows(use_pathway_list) %>%
        # 1. 创建 ONTOLOGY 因子，指定固定顺序
        mutate(ONTOLOGY = factor(ONTOLOGY, levels = c('BP', 'CC', 'MF', 'KEGG'))) %>%
        # 2. 按 ONTOLOGY 因子和 p.adjust 排序数据
        dplyr::arrange(ONTOLOGY, p.adjust) %>%
        # 3. 创建 Description 因子，决定Y轴顺序 (p值最小的在每个block的顶部)
        mutate(Description = factor(Description, levels = rev(unique(Description))))
      
      N_pathways <- nrow(use_pathway_group) # 总通路数，用于计算固定标签高度
      
      # 8. 准备左侧固定大小标签数据 (Requirement 1: BP->CC->MF->KEGG top-to-bottom, fixed size)
      print("    Preparing fixed-size labels...")
      label_categories <- factor(c('BP', 'CC', 'MF', 'KEGG'), levels = c('BP', 'CC', 'MF', 'KEGG'))
      total_categories <- length(label_categories)
      
      # 计算固定高度: 将总绘图高度 (约 N_pathways) 分成 N 等份
      # 每个标签占据 N_pathways / total_categories 的垂直空间
      block_size_y = ifelse(N_pathways > 0, N_pathways / total_categories, 1) 
      
      rect.data_group <- data.frame(ONTOLOGY = label_categories) %>%
        mutate(
          cat_index = as.numeric(ONTOLOGY), # BP=1, CC=2, MF=3, KEGG=4
          # 计算Y范围，留出固定比例的间隙 (e.g., 上下各留5% = 总共10%)
          gap_y = 0.1 * block_size_y, 
          # Y轴从下往上是 1 到 N_pathways
          # BP (index=1) 应该在最上面，对应Y轴数值大的区域
          # KEGG (index=4) 在最下面，对应Y轴数值小的区域
          # 所以需要反转 cat_index 的使用
          ymin_rect = (total_categories - cat_index) * block_size_y + 0.5 + gap_y, 
          ymax_rect = (total_categories - cat_index + 1) * block_size_y + 0.5 - gap_y,
          # 固定宽度
          xmin_rect = label_bar_x_start, 
          xmax_rect = label_bar_x_end
        ) #%>%
      # 不再过滤，总是显示所有4个标签
      # filter(ONTOLOGY %in% unique(use_pathway_group$ONTOLOGY)) 
      
      # 9. 绘制组合富集图 
      print("    Generating plot...")
      
      xaxis_max <- max(c(1, -log10(use_pathway_group$p.adjust)), na.rm = TRUE) + 1 
      
      p_combined <- use_pathway_group %>%
        ggplot(aes(x = -log10(p.adjust), y = Description, fill = ONTOLOGY)) + 
        # 条形图
        geom_col(width = 0.6, alpha = 0.8) + 
        # 通路描述文本
        geom_text(
          aes(x = 0.05, label = str_wrap(Description, width = 60)), 
          hjust = 0, size = 3, color = "black"
        ) +
        # 基因数量点图 (圆圈)
        geom_point(
          aes(x = gene_count_x_pos, size = Count), # 使用定义的X位置
          shape = 21, fill = "grey80", color = "black" 
        ) +
        # --- Requirement 3: Remove numbers inside circles ---
        # geom_text(
        #    aes(x = gene_count_x_pos, label = Count),
        #    size = 2.5, color = "black" 
        # ) +
        scale_size_continuous(name = 'Gene Count', range = c(2, 6)) + 
        # 左侧分类标签 (固定大小和顺序)
        geom_rect( 
          aes(xmin = xmin_rect, xmax = xmax_rect, ymin = ymin_rect, ymax = ymax_rect, fill = ONTOLOGY),
          data = rect.data_group,
          inherit.aes = FALSE,
          alpha = 0.8 
        ) +
        geom_text( # 标签文字
          aes(x = (xmin_rect + xmax_rect) / 2, y = (ymin_rect + ymax_rect) / 2, label = ONTOLOGY),
          data = rect.data_group,
          inherit.aes = FALSE,
          color = "black", size = 3.5 
        ) +
        # X轴轴线 
        geom_segment(
          aes(x = 0, y = 0.5, xend = xaxis_max, yend = 0.5), 
          linewidth = 1, inherit.aes = FALSE, color = "grey60"
        ) +
        # 标签和标度 
        labs(
          x = expression(-log[10]("Adjusted P-value")), 
          y = NULL, 
          title = paste0(dataset, " - ", group_name, ": Top GO & KEGG Pathways")
        ) +
        # 使用新的颜色方案 (Requirement 4)
        scale_fill_manual(name = 'Category', values = pal, 
                          # 图例顺序应与因子顺序一致 (BP在顶部)
                          guide = guide_legend(reverse = FALSE)) + 
        scale_x_continuous(
          breaks = scales::pretty_breaks(n = 4), 
          expand = expansion(mult = c(0.01, 0.05)) 
        ) +
        scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) + 
        # 主题 
        theme_prism(base_size = 10) + 
        theme(
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), 
          axis.line.x = element_line(color = "black"), axis.line.y = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.position = "right", 
          plot.title = element_text(hjust = 0.5, size=12, face="bold"), 
          plot.margin = margin(t=10, r=10, b=10, l=abs(label_bar_x_start)*8, unit = "pt") # 动态调整左边距
        ) +
        coord_cartesian(clip = 'off') # 允许画到绘图区域外
      
      # 10. 保存图
      plot_filename_base <- file.path(output_dir, paste0(dataset, "-", group_name, "_Enrichment_Combined"))
      num_pathways <- nrow(use_pathway_group)
      plot_height <- max(4, num_pathways * 0.25 + 1.5) 
      plot_width <- 9 
      
      ggsave(paste0(plot_filename_base, ".png"), plot = p_combined, height = plot_height, width = plot_width, limitsize = FALSE, bg = "white", dpi = 300)
      ggsave(paste0(plot_filename_base, ".pdf"), plot = p_combined, height = plot_height, width = plot_width, limitsize = FALSE, bg = "white")
      print(paste("    Combined plot saved for", group_name))
      
    } else {
      print(paste("    No significant pathways found to plot for", group_name))
    }
    
  } # 结束 group_name 循环
  
} # 结束 dataset 循环

print("All datasets processed.")