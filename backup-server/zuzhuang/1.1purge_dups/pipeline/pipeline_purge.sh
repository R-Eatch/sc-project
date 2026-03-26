#!/bin/bash
export PATH=/public/home/zhengjiangmin/Software/minimap2-2.22_x64-linux:$PATH
export PATH=/public/home/sunxuebo/software/purge_dups-master/bin:$PATH
# 设置主要组装的路径
pri_asm="../genome.nextpolish.fasta"

# 设置PacBio数据文件列表的路径
pb_data="/data01/sunxuebo/project/zuzhuang/1.hifiasm/hifi_reads.fastq"
# 第1步：对于PacBio CLR读数


minimap2 -I 4G  -xmap-hifi $pri_asm "$pb_data" | gzip -c - >  "hifi_reads.paf.gz"

# 生成PB.base.cov和PB.stat文件
pbcstat *.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log

# 注意：如果您有大型基因组，请设置minimap2 -I选项

# 第1步：分割组装并进行自我对齐
split_fa $pri_asm > $pri_asm.split
minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

# 第2步：清除haplotigs和重叠
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

# 第3步：从草稿组装中获取清洁的主要和haplotig序列
get_seqs -e dups.bed $pri_asm

# 第4步：合并hap.fa和$hap_asm，重做以上步骤以获得合适的haplotig集
# （如果您有备用组装）

# 限制：对重复序列的处理能力有限
