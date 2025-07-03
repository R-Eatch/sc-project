# ① 安装（一次即可）
install.packages("babelgene")   # <1 MB，纯 CRAN，离线数据内置

# ② 加载
library(babelgene)
library(data.table)

# ③ 读入你的 CSV
csv_in  <- "H-MG_ranked_genes_author_cell_type.csv"
csv_out <- "H-MG_ranked_genes_author_cell_type-mouse.csv"
df <- fread(csv_in)

# ④ 批量查询：ENS* → 小鼠基因名
#    orthologs() 自动识别 Ensembl/Entrez/符号；human=TRUE 表示输入基因属人类
mouse_map <- orthologs(
  genes   = unique(df$names),   # 提供所有待转换的 ID
  species = "mouse",            # 输出物种
  human   = TRUE                # 输入物种（人）
)

# ⑤ 构建映射字典（人 ENSG → 小鼠 symbol）
dict <- setNames(mouse_map$symbol, mouse_map$human_ensembl)

# ⑥ 替换；无同源则保留原 ID
df[, mouse_gene := fifelse(names %in% names(dict), dict[names], names)]

# ⑦ 保存
fwrite(df, csv_out)
cat("✅ babelgene 转换完成：", csv_out, "\n")
