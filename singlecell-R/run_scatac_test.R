library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(rtracklayer)
counts <- Read10X_h5(filename = "D:/BaiduNetdiskDownload/R-LA-MG-ATAC/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "D:/BaiduNetdiskDownload/R-LA-MG-ATAC/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "D:/BaiduNetdiskDownload/R-LA-MG-ATAC/fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)

pbmc[['peaks']]
granges(pbmc)
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc <- pbmc[as.vector(peaks.keep), ]



# --- 1. 定义输入和参数 ---
gtf_file <- "D:/111/motif/genome/rabbit.gtf"  # *** 你的 GTF 文件路径 ***
genome_version <- "orycun2"                   # *** 你的基因组版本标识符 ***
feature_type_to_use <- "transcript"          # *** 推荐使用 "transcript" ***

# (可选) 如果你的 GTF 有生物类型信息，并且你想筛选，在这里指定
# biotypes_to_keep <- c("protein_coding", "lincRNA", "processed_transcript") 
# gtf_biotype_attribute <- "transcript_biotype" # GTF第9列中代表生物类型的属性名 (需要你自己确认！)

# 假设 pbmc 对象已存在
if (!exists("pbmc") || !inherits(pbmc, "Seurat")) {
  stop("错误：未能找到名为 'pbmc' 的 Seurat 对象。")
}

# --- 2. 导入 GTF 文件 ---
message("正在导入 GTF 文件: ", gtf_file)
tryCatch({
  gtf_all <- rtracklayer::import(gtf_file, format = "gtf")
}, error = function(e) {
  stop("导入 GTF 文件失败: ", e$message)
})
message("GTF 导入成功，总特征数: ", length(gtf_all))

# --- 3. 筛选主要特征类型 (推荐 transcript) ---
message("正在根据类型 '", feature_type_to_use, "' 筛选特征...")
annotation_gr <- gtf_all[mcols(gtf_all)$type == feature_type_to_use]

if (length(annotation_gr) == 0) {
  stop("根据类型 '", feature_type_to_use, "' 筛选后得到 0 个特征。请检查 GTF 文件或尝试其他类型 (如 'gene')。")
}
message("筛选后得到 ", length(annotation_gr), " 个类型为 '", feature_type_to_use, "' 的特征。")

# --- 4. 筛选标准染色体 ---
message("正在筛选标准染色体...")
# 获取 pbmc 对象中 peaks 的标准染色体 （作为参考）
# 注意：standardChromosomes() 可能对你的基因组标识符 'orycun2' 无效，
# 如果无效，我们需要手动定义或直接使用 pbmc 对象中的染色体
pbmc_seqlevels <- seqlevels(granges(pbmc))
standard_chrs <- tryCatch(standardChromosomes(annotation_gr), error = function(e) NULL)

if (!is.null(standard_chrs)) {
  message("使用 standardChromosomes() 识别标准染色体。")
  chrs_to_keep <- intersect(standard_chrs, seqlevels(annotation_gr))
} else {
  message("standardChromosomes() 可能不支持此基因组。将使用 pbmc 对象中的染色体作为标准。")
  # 确保只保留在 pbmc 对象中也存在的染色体
  chrs_to_keep <- intersect(pbmc_seqlevels, seqlevels(annotation_gr))
}

if (length(chrs_to_keep) == 0) {
  warning("未能确定要保留的标准染色体，或 pbmc 与 GTF 染色体无交集。跳过染色体筛选。请检查两者染色体名称。")
} else {
  message("将保留以下染色体: ", paste(head(chrs_to_keep), collapse=", "), "...")
  original_len <- length(annotation_gr)
  seqlevels(annotation_gr, pruning.mode="coarse") <- chrs_to_keep # 保留指定染色体上的特征
  # 某些情况下，如果GRanges中包含不在chrs_to_keep列表中的序列，需要先过滤
  # annotation_gr <- annotation_gr[seqnames(annotation_gr) %in% chrs_to_keep]
  # seqlevels(annotation_gr) <- chrs_to_keep # 更新seqlevels列表
  message("筛选掉 ", original_len - length(annotation_gr), " 个位于非标准或非 pbmc 染色体上的特征。")
}


# --- 5. (可选) 筛选生物类型 ---
# 只有当你确定 GTF 文件中有可靠的生物类型信息时才执行此步
# if (exists("biotypes_to_keep") && exists("gtf_biotype_attribute")) {
#   message("正在根据生物类型筛选...")
#   mcols_data <- mcols(annotation_gr)
#   if (gtf_biotype_attribute %in% colnames(mcols_data)) {
#     original_len <- length(annotation_gr)
#     keep_idx <- mcols_data[[gtf_biotype_attribute]] %in% biotypes_to_keep
#     # 处理可能的 NA 值
#     keep_idx[is.na(keep_idx)] <- FALSE
#     annotation_gr <- annotation_gr[keep_idx]
#     message("根据生物类型筛选后，保留 ", length(annotation_gr), " 个特征 (移除了 ", original_len - length(annotation_gr), ")。")
#   } else {
#     warning("指定的生物类型属性 '", gtf_biotype_attribute, "' 在 GTF 元数据中未找到，跳过按生物类型筛选。")
#   }
# } else {
#   message("未指定生物类型筛选条件，将保留所有类型的 '", feature_type_to_use, "' 特征。")
# }

# --- 6. 确保 'gene_name' 元数据列存在且尽可能填充 ---
message("正在处理 'gene_name' 元数据列...")
mcols_data <- mcols(annotation_gr)
available_cols <- colnames(mcols_data)

# 确定填充 NA 的备用 ID 列 (优先 gene_id, 其次 transcript_id)
id_to_use <- NA
if ("gene_id" %in% available_cols) {
  id_to_use <- "gene_id"
  message("将优先使用 'gene_id' 填充缺失的 'gene_name'。")
} else if ("transcript_id" %in% available_cols && feature_type_to_use == "transcript") {
  id_to_use <- "transcript_id" # 只有当我们在处理转录本时，用 transcript_id 才有意义
  message("未找到 'gene_id'，将使用 'transcript_id' 填充 'gene_name'。")
} else {
  warning("未找到 'gene_id' 或合适的备用 ID ('transcript_id') 来填充 'gene_name'。")
}

# 确保 'gene_name' 列存在
if (!"gene_name" %in% available_cols) {
  message("'gene_name' 列不存在，将创建此列。")
  mcols(annotation_gr)$gene_name <- NA_character_
  mcols_data <- mcols(annotation_gr) # 更新视图
}

# 填充 NA 值
if (!is.na(id_to_use)) {
  indices_to_fix <- which(is.na(mcols_data$gene_name) & !is.na(mcols_data[[id_to_use]]))
  num_to_fix <- length(indices_to_fix)
  if (num_to_fix > 0) {
    message("发现 ", num_to_fix, " 个特征的 'gene_name' 为 NA。正在使用 '", id_to_use, "' 进行填充...")
    mcols(annotation_gr)$gene_name[indices_to_fix] <- mcols_data[[id_to_use]][indices_to_fix]
    message("填充完成。")
  } else {
    message("没有找到需要使用 '", id_to_use, "' 修复的 'gene_name'。")
  }
}

# 最终检查 'gene_name'
final_non_na_gene_names <- sum(!is.na(mcols(annotation_gr)$gene_name))
if (length(annotation_gr) > 0 && final_non_na_gene_names == 0) {
  warning("警告：处理后所有特征的 'gene_name' 仍为 NA！下游分析可能受影响。请检查 GTF 文件属性和脚本逻辑。")
} else {
  message("最终具有非 NA 'gene_name' 的特征数量: ", final_non_na_gene_names, " / ", length(annotation_gr))
}
# 检查重复 gene_name
num_duplicates <- sum(duplicated(mcols(annotation_gr)$gene_name[!is.na(mcols(annotation_gr)$gene_name)]))
if(num_duplicates > 0) {
  warning(num_duplicates, " 个重复的 'gene_name' 值存在（在非 NA 值中）。")
}

# --- 7. 统一染色体命名风格并设置基因组信息 ---
message("正在统一染色体命名风格并设置基因组信息...")
# 获取 pbmc 对象的风格
pbmc_style <- tryCatch(seqlevelsStyle(granges(pbmc))[1], error = function(e) NA)

if (!is.na(pbmc_style)) {
  message("检测到 pbmc 对象的风格为: ", pbmc_style, ". 将尝试统一注释对象的风格。")
  # 尝试应用风格 (这会自动处理命名转换，如 1 <-> chr1)
  # 需要在分配 genome() 之前或同时进行
  tryCatch({
    seqlevelsStyle(annotation_gr) <- pbmc_style
    message("注释对象的染色体命名风格已设置为: ", pbmc_style)
  }, error = function(e) {
    warning("尝试自动统一染色体命名风格失败: ", e$message, " 请确保手动检查两者风格一致。")
  })
} else {
  warning("无法确定 pbmc 对象的染色体命名风格。请务必手动确认注释对象与 pbmc 对象的染色体命名一致！")
}

# 设置基因组版本
tryCatch({
  genome(annotation_gr) <- genome_version
  message("注释对象的基因组已设置为: ", genome_version)
}, error = function(e) {
  warning("设置基因组信息时出错: ", e$message)
})

# --- 8. 最终检查并分配注释 ---
message("最终检查处理后的注释对象:")
print(annotation_gr)
print("元数据列:")
print(colnames(mcols(annotation_gr)))

message("将处理后的注释分配给 Seurat 对象 'pbmc'...")
tryCatch({
  Annotation(pbmc) <- annotation_gr
  message("注释成功分配给 'pbmc'。")
  message("已分配注释的前几行:")
  print(head(Annotation(pbmc)))
}, error = function(e) {
  stop("将注释分配给 Seurat 对象时发生错误: ", e$message)
})

message("注释处理和分配流程完成。")


# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc)
DensityScatter(pbmc, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# add fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
#FragmentHistogram(object = pbmc)
VlnPlot(
  object = ea,
  features = c('nCount_peaks', 'TSS.enrichment', 'nucleosome_signal', 'pct_reads_in_peaks'),
  pt.size = 0.1,
  ncol = 4
)
ea_raw <-pbmc

ea <- subset(
  x = ea_raw,
  subset = nCount_peaks > 9000 &
    nCount_peaks < 100000 &
    pct_reads_in_peaks > 40 &
    #blacklist_ratio < 0.01 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)
saveRDS(ea,'R-MG-LA-normal.rds')
ea <- RunTFIDF(ea)
ea <- FindTopFeatures(ea, min.cutoff = 'q0')
ea <- RunSVD(ea)
ea <- RunUMAP(object = ea, reduction = 'lsi', dims = 2:30)
ea <- FindNeighbors(object = ea, reduction = 'lsi', dims = 2:30)
ea <- FindClusters(object = ea, verbose = FALSE, algorithm = 3)
DimPlot(object = ea, label = TRUE) + NoLegend()
gene.activities <- GeneActivity(ea)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
ea[['RNA']] <- CreateAssayObject(counts = gene.activities)
ea <- NormalizeData(
  object = ea,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(ea$nCount_RNA)
)
DefaultAssay(ea) <- 'RNA'
ea <- ScaleData(ea)
ea <- RunPCA(ea)
ea <- RunUMAP(object = ea, reduction = 'pca', dims = 1:30)
ea <- FindNeighbors(object = ea, reduction = 'pca', dims = 1:30)
ea <- FindClusters(object = ea, verbose = FALSE, algorithm = 3)
DimPlot(object = ea, label = TRUE) + NoLegend()
FeaturePlot(
  object = ea,
  features = c('LYZG', 'ESR1', 'PRLR', 'ACTA2', 'CSN2', 'CSN3'),
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
DefaultAssay(ea) <- 'peaks'
Idents(ea)
ea <- SortIdents(ea)

CoveragePlot(
  object = ea,
  region = "LYZ",
  extend.upstream = 1000,
  extend.downstream = 1000
)

