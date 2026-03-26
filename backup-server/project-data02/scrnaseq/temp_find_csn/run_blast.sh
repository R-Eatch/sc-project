#!/usr/bin/env bash
# ------------------------------------------------------------
#  find_csn_in_sugar_glider.sh
#  自动比对 CSN 蛋白到蜜袋鼯基因组，并交叉注释
#  用法:
#    bash find_csn_in_sugar_glider.sh  csn.fa  genome.fa  genome.gtf  [outprefix]
# ------------------------------------------------------------
set -euo pipefail

###########################
# 参数检查
###########################
if [[ $# -lt 3 ]]; then
  echo "用法: $0  <csn.fa>  <genome.fa>  <genome.gtf>  [outprefix]" >&2
  exit 1
fi

CSN_FASTA="$(realpath "$1")"
GENOME_FA="$(realpath "$2")"
GTF_FILE="$(realpath "$3")"
OUTPREFIX="${4:-csn_vs_glider}"

###########################
# 1. 建立 BLAST 数据库
###########################
echo "▶ 建立 BLAST 数据库..."
if [[ ! -f "${GENOME_FA}.nin" ]]; then      # 已存在则跳过
  makeblastdb -in "$GENOME_FA" -dbtype nucl -out "${GENOME_FA}"
fi

###########################
# 2. 运行 tblastn
###########################
echo "▶ 运行 tblastn..."
tblastn \
  -query  "$CSN_FASTA" \
  -db     "$GENOME_FA" \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -num_threads "$(nproc --ignore=20)" \
  -out    "${OUTPREFIX}.tblastn.tsv"

###########################
# 3. 将 tblastn 结果转成 BED
###########################
echo "▶ 转换为 BED 格式..."
awk 'BEGIN{OFS="\t"}{
      start = ($9 < $10 ? $9 : $10);         # sstart/send 取较小者
      end   = ($9 > $10 ? $9 : $10);
      print $2, start-1, end, $1, $13, ($9<$10?"+":"-")
    }' "${OUTPREFIX}.tblastn.tsv" \
  | sort -k1,1 -k2,2n > "${OUTPREFIX}.hits.bed"

###########################
# 4. 与 GTF 注释交叉
###########################
echo "▶ 与 GTF 交叉..."
bedtools intersect -a "${OUTPREFIX}.hits.bed" -b "$GTF_FILE" -wa -wb \
  > "${OUTPREFIX}.gtf_overlap.tsv"

echo "▶ 提取重叠 gene_id 列表..."
grep -o 'gene_id "[^"]*"' "${OUTPREFIX}.gtf_overlap.tsv" \
  | awk '{print $2}' | tr -d '"' \
  | sort -u > "${OUTPREFIX}.putative_geneIDs.txt"

###########################
# 5. 完成
###########################
echo "✅ 运行完毕！
  - BLAST 结果:        ${OUTPREFIX}.tblastn.tsv
  - 命中 BED:          ${OUTPREFIX}.hits.bed
  - GTF 重叠注释:      ${OUTPREFIX}.gtf_overlap.tsv
  - 候选 CSN 基因 ID:  ${OUTPREFIX}.putative_geneIDs.txt"

