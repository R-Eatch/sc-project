bedtools intersect -a output/filtered_alignment_fixed.bed -b /data01/sunxuebo/project/atac/2.liftover/J103660_S-GES14-MG-bATAC_peaks.narrowPeak -wa -u -f 0.5 > intersected_peaks.bed
##建议手动复制然后运行，不知道为什么bedtools无法调用
