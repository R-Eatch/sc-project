#!/usr/bin/env bash
###############################################################################
#  map_csn_rnaseq.sh
#  从 CSN tblastn 命中窗口提取转录本 → 生成补充注释
#
#  依赖：bedtools, minimap2, samtools, gffread
###############################################################################

set -euo pipefail

###################### ① 路径配置 ############################################
# ——按需修改——
GENOME="/data01/sunxuebo/project/zuzhuang/final_genome_gff/sugar_glider.fa"
TRANSCRIPT_FA="/data01/sunxuebo/project/zuzhuang/00.translig/trans_out/MG/TransLiG_Out_Dir/TransLiG.fa"
HITS_BED="/data02/sunxuebo/project/scrnaseq/temp_find_csn/sg.hits.bed"           # tblastn 命中的原始 BED
WORKDIR="csn_pipeline_out"       # 输出目录
THREADS=16                       # 并行线程数
SLOP=15000                       # tblastn 区域外扩 (bp)

mkdir -p "${WORKDIR}"
cd       "${WORKDIR}"

###################### ② 取命中窗口 ±15 kb ###################################
echo "▶ Step 1: slop BED..."
bedtools slop \
  -i  "${HITS_BED}" \
  -g  "${GENOME}.fai" \
  -b  "${SLOP}" > csn_hits_slop.bed

###################### ③ minimap2 转录本比对并截取窗口 #######################
echo "▶ Step 2: minimap2 mapping..."
minimap2 -ax splice -uf --secondary=no -t "$THREADS" "$GENOME" "$TRANSCRIPT_FA" \
  | samtools sort -o rnaseq_splice.coord.bam
samtools index rnaseq_splice.coord.bam

# (2) intersect 只截窗口
bedtools intersect -a rnaseq_splice.coord.bam \
                   -b csn_hits_slop.bed -wa -ubam \
  > csn_hits_transcripts.bam

# (3) 名称排序 → gffread
samtools sort -n -o csn_hits_transcripts.name.bam csn_hits_transcripts.bam
gffread csn_hits_transcripts.name.bam -g "$GENOME" \
        -T -o csn_from_rnaseq.gtf -w csn_from_rnaseq.fa
###################### ⑤ (可选) 坐标排序 & index #############################
samtools sort -o csn_hits_transcripts.coord.bam csn_hits_transcripts.bam
samtools index csn_hits_transcripts.coord.bam

echo "✅ 完成：
  • 比对 BAM：           ${WORKDIR}/csn_hits_transcripts.coord.bam
  • 补充注释 GTF：       ${WORKDIR}/csn_from_rnaseq.gtf
  • 转录本 FASTA：       ${WORKDIR}/csn_from_rnaseq.fa"
###############################################################################

