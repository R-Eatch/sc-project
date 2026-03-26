#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
from matplotlib.sankey import Sankey

###############################################################################
# 1) 在这里硬编码三个文件的路径
###############################################################################
GFF_FILE = "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20240230.gff"       # 人类调控元件 GFF
BED1_FILE = "./output/filtered_alignment_fixed.bed"  # 比对到蜜袋鼯基因组后的 BED
BED2_FILE = "./intersected_peaks.bed"                # 与 ATAC-seq 交集后的 BED

###############################################################################
# 2) 解析函数：统计各文件中的调控元件类型及数量
###############################################################################
def parse_gff(gff_file):
    """
    解析 GFF 文件，统计第三列（feature_type）中各种调控元件出现次数。
    返回字典 { 'enhancer': 数量, 'promoter': 数量, ... }
    """
    counts = {}
    with open(gff_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            # 基本格式检查
            if len(fields) < 3:
                continue
            feature_type = fields[2]  # GFF 第三列通常是 feature，如 "enhancer", "promoter" 等
            counts[feature_type] = counts.get(feature_type, 0) + 1
    return counts

def parse_bed(bed_file):
    """
    解析 BED 文件，假设第 4 列包含类似 "enhancer:xxx" 或 "promoter:xxx" 等信息。
    若存在冒号，则取冒号前的部分作为类型；否则直接用整列作类型。
    返回字典 { 'enhancer': 数量, 'promoter': 数量, ... }
    """
    counts = {}
    with open(bed_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            annotation = fields[3]
            # 取 ":" 前的部分作为类型，否则用完整字段
            if ':' in annotation:
                feature_type = annotation.split(':')[0]
            else:
                feature_type = annotation

            counts[feature_type] = counts.get(feature_type, 0) + 1
    return counts

###############################################################################
# 3) 柱状图绘制函数（柱子上方显示数值，并输出到 ./output）
###############################################################################
def plot_bar_chart(counts, title, output_prefix):
    """
    绘制柱状图，每个柱子顶部显示数值。
    counts: dict, {类型: 数量}
    title: 图的标题（英文）
    output_prefix: 输出文件名前缀
    """
    labels = list(counts.keys())
    values = [counts[l] for l in labels]

    plt.figure()
    bars = plt.bar(labels, values)
    plt.title(title)
    plt.xlabel("Regulatory Element Type")
    plt.ylabel("Count")
    plt.tight_layout()

    # 在柱顶显示数值
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, height,
                 f"{int(height)}",
                 ha='center', va='bottom')

    # 保存图像到 ./output
    output_file = f"./output/{output_prefix}_bar_chart.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"[INFO] Bar chart saved to: {output_file}")

###############################################################################
# 4) Sankey 图绘制函数
#    设定：假设 BED1 是 GFF 的子集，BED2 是 BED1 的子集
###############################################################################
def plot_sankey(gff_total, bed1_total, bed2_total):
    """
    绘制 Sankey 图，展示三个数据集之间的层级关系:
      - GFF 全部
        -> 一部分映射到 BED1, 另一部分未映射
      - BED1
        -> 一部分与 ATAC 交集 (BED2), 另一部分不在交集
    若 bed1_total > gff_total / bed2_total > bed1_total，会导致流向不符。
    """
    fig, ax = plt.subplots()
    sankey = Sankey(ax=ax, unit="", format="%.0f", headangle=120, scale=1.0, gap=0.5)

    # 第一个节点：GFF => BED1 + 未映射
    # flows = [ +GFF, -BED1, -(GFF - BED1)]
    sankey.add(
        flows=[gff_total, -bed1_total, -(gff_total - bed1_total)],
        labels=["GFF", "BED1", "Not mapped"],
        orientations=[0, 0, 0],
    )

    # 第二个节点：BED1 => BED2 + 不在交集
    # flows = [ +BED1, -BED2, -(BED1 - BED2)]
    # 与上一个节点相连 (prior=0), 连接第1个节点的 flow index=1 (即 BED1)
    sankey.add(
        flows=[bed1_total, -bed2_total, -(bed1_total - bed2_total)],
        labels=[None, "BED2", "Not in intersection"],
        orientations=[0, 0, 0],
        prior=0,       # 与上一个节点索引 0 相连
        connect=(1, 0) # 本节点的第0号流 -> 上一个节点的第1号流
    )

    sankey_diagrams = sankey.finish()
    ax.set_title("Sankey Diagram of Three Datasets")
    output_file = "./output/sankey_diagram.png"
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"[INFO] Sankey diagram saved to: {output_file}")

###############################################################################
# 5) 主函数：执行统计、输出信息、绘图
###############################################################################
def main():
    # 1) 统计人类 GFF
    gff_counts = parse_gff(GFF_FILE)
    gff_total = sum(gff_counts.values())
    print(f"\n[ GFF ]: {GFF_FILE}")
    print(f"  Total regulatory elements: {gff_total}")
    for k, v in gff_counts.items():
        print(f"  {k}: {v}")

    # 2) 统计比对到蜜袋鼯基因组后的 BED
    bed1_counts = parse_bed(BED1_FILE)
    bed1_total = sum(bed1_counts.values())
    print(f"\n[ BED1 ]: {BED1_FILE}")
    print(f"  Total regulatory elements: {bed1_total}")
    for k, v in bed1_counts.items():
        print(f"  {k}: {v}")

    # 3) 统计与 ATAC-seq 交集后的 BED
    bed2_counts = parse_bed(BED2_FILE)
    bed2_total = sum(bed2_counts.values())
    print(f"\n[ BED2 ]: {BED2_FILE}")
    print(f"  Total regulatory elements: {bed2_total}")
    for k, v in bed2_counts.items():
        print(f"  {k}: {v}")

    # 绘制柱状图
    plot_bar_chart(gff_counts,
                   "Human Regulatory Elements Statistics",
                   "gff")
    plot_bar_chart(bed1_counts,
                   "Mapped to Sugar Glider Genome - Regulatory Elements Statistics",
                   "bed1")
    plot_bar_chart(bed2_counts,
                   "Specified Sample Intersection - Regulatory Elements Statistics",
                   "bed2")

    # 绘制 Sankey 图（若数据符合子集关系）
    plot_sankey(gff_total, bed1_total, bed2_total)

if __name__ == "__main__":
    main()

