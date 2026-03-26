#!/bin/bash

# --- 配置 ---
# 遇到任何错误则停止脚本
set -euo pipefail

# 输入文件
READS_FASTQ="/data01/sunxuebo/project/zuzhuang/1.hifiasm/hifi_reads.fastq"  # !!! 修改为您的 FASTQ 文件名 (可以是 .fastq 或 .fastq.gz)
GENOME_FASTA="/data01/sunxuebo/project/zuzhuang/final_genome_gff/sugar_glider.fa" # !!! 修改为您的基因组文件名

# 输出文件 (中间和最终结果)
SORTED_BAM_FILE="alignment.sorted.bam" # 比对、排序后的BAM文件名
STATS_FILE="mapping_stats.txt"       # 存放统计结果的文件名

# 使用的线程数 (根据您的服务器核心数调整)
THREADS=8

# minimap2 比对预设参数 (重要!)
# - 'map-hifi': 推荐用于 PacBio HiFi reads (即使已转为 FASTQ)
# - 'map-pb': 推荐用于 PacBio CLR reads (错误率较高)
# - 'map-ont': 推荐用于 Oxford Nanopore reads
# - 'sr': 推荐用于短读长 (Illumina) Paired-End reads
# !!! 根据您的 reads 类型选择合适的预设 !!!
MINIMAP2_PRESET="map-hifi"  # 假设您的 reads 源自 HiFi

# --- 检查 ---
echo "--- 开始计算映射率（挂载率） Pipeline (FASTQ 输入) ---"

# 检查 minimap2 和 samtools 是否安装
if ! command -v minimap2 &> /dev/null; then
    echo "错误: 未找到 minimap2 命令。请先安装 minimap2。"
    exit 1
fi
if ! command -v samtools &> /dev/null; then
    echo "错误: 未找到 samtools 命令。请先安装 samtools。"
    exit 1
fi
echo "[信息] minimap2 已找到: $(command -v minimap2)"
echo "[信息] samtools 已找到: $(command -v samtools)"

# 检查输入文件是否存在
if [ ! -f "$READS_FASTQ" ]; then
    echo "错误: 输入的 FASTQ 文件 '$READS_FASTQ' 不存在！"
    exit 1
fi
if [ ! -f "$GENOME_FASTA" ]; then
    echo "错误: 输入的基因组文件 '$GENOME_FASTA' 不存在！"
    exit 1
fi
echo "[信息] 使用的 Reads 文件: $READS_FASTQ"
echo "[信息] 使用的参考基因组: $GENOME_FASTA"
echo "[信息] 使用 minimap2 预设参数: $MINIMAP2_PRESET"
echo "[信息] 使用线程数: $THREADS"

# --- Pipeline 步骤 ---

# 1. 比对 reads 到基因组 (minimap2)，并将 SAM 输出通过管道传给 samtools 进行 BAM 转换和排序
echo "[步骤 1 & 2] 正在使用 minimap2 比对并将结果转换为排序后的 BAM 文件..."
minimap2 -ax "$MINIMAP2_PRESET" -t $THREADS "$GENOME_FASTA" "$READS_FASTQ" | samtools view -@ $THREADS -bS - | samtools sort -@ $THREADS - -o "$SORTED_BAM_FILE"
#  - minimap2: 进行比对, 输出 SAM 到标准输出 (stdout)
#  - samtools view -@ $THREADS -bS - : 从标准输入 (stdin, 由 '-' 指定) 读取 SAM, 转换为 BAM (-bS), 使用多线程, 输出 BAM 到 stdout
#  - samtools sort -@ $THREADS - -o "$SORTED_BAM_FILE" : 从 stdin ('-') 读取 BAM, 进行排序, 使用多线程, 将排序后的 BAM 写入输出文件 (-o)
echo "[步骤 1 & 2] 比对、转换和排序完成。排序后的 BAM 文件: $SORTED_BAM_FILE"

# 3. 为排序后的 BAM 文件创建索引
echo "[步骤 3] 正在为排序后的 BAM 文件创建索引..."
samtools index -@ $THREADS "$SORTED_BAM_FILE"
# 索引文件会自动命名为 ${SORTED_BAM_FILE}.bai
echo "[步骤 3] BAM 文件索引完成。"

# 4. 使用 samtools flagstat 计算比对统计信息
echo "[步骤 4] 正在使用 samtools flagstat 计算比对统计信息..."
samtools flagstat -@ $THREADS "$SORTED_BAM_FILE" > "$STATS_FILE"
echo "[步骤 4] 统计信息计算完成。结果已保存至 '$STATS_FILE'。"

# 5. 从统计文件中提取并打印映射率
echo "[步骤 5] 正在从 '$STATS_FILE' 文件中提取映射率..."
# 使用 grep 找到包含 "mapped" 的行，然后用 awk 提取括号中的百分比
MAPPED_LINE=$(grep -E '^[0-9]+ \+ [0-9]+ mapped \(' "$STATS_FILE")
# 使用 awk，以 '(' 或 '%' 或 ':' 作为分隔符，提取第二个字段
MAPPING_RATE=$(echo "$MAPPED_LINE" | awk -F'[():%]' '{print $2}')

if [ -n "$MAPPING_RATE" ]; then
    echo "--------------------------------------------------"
    echo " 映射率 (挂载率): ${MAPPING_RATE}%"
    echo "--------------------------------------------------"
    echo "排序后的比对文件: $SORTED_BAM_FILE"
    echo "完整的统计信息请查看文件: $STATS_FILE"
else
    echo "错误: 未能自动从 '$STATS_FILE' 文件中提取映射率。"
    echo "请手动检查该文件。"
fi

echo "--- Pipeline 执行完毕 ---"

exit 0
