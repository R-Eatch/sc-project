library(ArchR)
set.seed(1)
addArchRThreads(threads = 1) 
addArchRChrPrefix(chrPrefix = FALSE)
library(rtracklayer)
library(GenomicRanges)
library(dplyr)


# 读取 GTF 文件
gtf <- import("D:/111/motif/genome/rabbit.gtf")  # 替换为你的路径
gtf <- gtf[gtf$type %in% c("gene", "exon"), ]

# 提取 gene
genes <- gtf[gtf$type == "gene"]
mcols(genes)$symbol <- genes$gene_name  # 确保有名为 symbol 的列

# 提取 exons
exons <- gtf[gtf$type == "exon"]
mcols(exons)$symbol <- exons$gene_name  # 同样要有 symbol 列

# 提取 TSS（转录起始位点）
tss <- genes
start(tss) <- ifelse(strand(tss) == "+", start(tss), end(tss))
end(tss) <- start(tss)
tss$symbol <- genes$symbol
geneAnnotation <- createGeneAnnotation(
  genes = genes,
  exons = exons,
  TSS = tss
)
library(Biostrings)
library(GenomicRanges)

# 读取你的 genome.fa 文件
genome_fa <- readDNAStringSet("D:/111/motif/genome/rabbit_CLEAN.fa")  # 替换路径


# 提取染色体名称和长度
chromSizes <- seqlengths(genome_fa)
# 创建一个空的 GRanges 对象作为 blacklist（可选，但推荐）
blacklist <- GRanges()
genomeAnnotation <- list(
  genome = genome_fa,
  chromSizes = chromSizes,
  blacklist = blacklist
)
inputFiles <- c("D:/BaiduNetdiskDownload/R-LA-MG-ATAC/fragments.tsv.gz")
sampleNames <- c("R-MG-LA")  # 与 inputFiles 一一对应

# 可选：指定输出路径（默认当前目录）
outputDir <- "D:/111/"
keep_chr <- c(as.character(1:21), "X", "MT")   # 与 clean_fa 一致

geneAnnotation_fix <- list(
  genes = keepSeqlevels(geneAnnotation$genes, keep_chr, pruning.mode = "coarse"),
  exons = keepSeqlevels(geneAnnotation$exons, keep_chr, pruning.mode = "coarse"),
  TSS   = keepSeqlevels(geneAnnotation$TSS,   keep_chr, pruning.mode = "coarse")
)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = 'sample1',
  geneAnnotation = geneAnnotation_fix,             # 你创建的对象
  genomeAnnotation = genomeAnnotation,
  logFile ="D:/111/arrow_debug.log",
  force = TRUE# 你创建的对象
)


