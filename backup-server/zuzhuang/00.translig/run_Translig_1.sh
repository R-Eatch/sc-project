#!/bin/bash

# 样本来源目录
sample_source_dir="/data01/sunxuebo/project/zuzhuang/transcriptome_data/ANNO_XS01KF2023120164_PM-XS01KF2023120164-03_2024-01-02_15-14-29_222M5HLT4/Cleandata"

# 样本结果输出目录
output_dir="/data01/sunxuebo/project/zuzhuang/00.translig/trans_out"

# 组装命令的路径
translig_cmd="/public/home/sunxuebo/software/TransLiG_1.3/TransLiG"

# 遍历样本来源目录中的每个子目录
for sample_folder in $sample_source_dir/*; do
    if [ -d "$sample_folder" ]; then
        sample_name=$(basename "$sample_folder")

        # 创建PBS任务文件
        pbs_file="translig_${sample_name}.pbs"
        cat << EOF > $pbs_file
#PBS -N translig_$sample_name
#PBS -o /data01/sunxuebo/project/zuzhuang/00.translig/
#PBS -e /data01/sunxuebo/project/zuzhuang/00.translig/
#PBS -q more
#PBS -l nodes=1:ppn=20
#PBS -r y

mkdir -p $output_dir/$sample_name
cd $output_dir/$sample_name
r1_file="${sample_source_dir}/${sample_name}/${sample_name}_R1.fq.gz"
r2_file="${sample_source_dir}/${sample_name}/${sample_name}_R2.fq.gz"

$translig_cmd -s fq -p pair -l \$r1_file -r \$r2_file -m RF -k 31 -g 200
EOF

        # 提交PBS任务
        qsub $pbs_file
    fi
done

