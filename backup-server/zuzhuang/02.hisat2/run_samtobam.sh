#!/bin/bash

# 创建输出目录，如果不存在的话
mkdir -p bamfile
mkdir -p finalbam
# 转换SAM文件为BAM文件
for dir in hisat2_out/*; do
    if [ -d "$dir" ]; then
        organ_name=$(basename "$dir")
        echo "Processing $organ_name"
        samtools view -@ 32 -bS "$dir/$organ_name.sam" > "bamfile/output_$organ_name.bam"
    fi
done

# 合并BAM文件
samtools merge -@ 32 finalbam/merged.bam bamfile/output_*.bam

echo "Merge completed"


