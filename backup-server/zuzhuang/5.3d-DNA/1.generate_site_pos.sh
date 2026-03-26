#!/usr/bin/env bash

# 1. generate site.pos

bwa="/public/home/zhengjiangmin/Software/bwa/bwa"
samtools="/public/home/zhengjiangmin/Software/samtools/bin/samtools"
bioawk="/public/home/zhengjiangmin/Software/bioawk-1.0/bioawk"
juicer_path="/public/home/zhengjiangmin/Software/juicer-1.6"
enzyme_type="MboI"
export PATH=/public/home/zhengjiangmin/Software/bwa:$PATH



rawGenome="sugar_glider.fa"



[[ -f ${rawGenome}.pac ]] || $bwa index ${rawGenome}
[[ -f ${rawGenome}.fai ]] || $samtools faidx ${rawGenome}

# field1: restriction enzyme type
# field2: genome name
# field3: genome path
[[ -f `basename $rawGenome`_$enzyme_type.txt ]] \
 || python3 ${juicer_path}/misc/generate_site_positions.py $enzyme_type `basename ${rawGenome}` ${rawGenome}

# 2. generate genome size file
[[ -f ${rawGenome}.sizes ]] \
 || $bioawk -c fastx '{print $name, length($seq)}' ${rawGenome} > ${rawGenome}.sizes

# 3. RUN juicer
/public/home/sunxuebo/software/juicer-1.6/CPU/juicer.sh -g `basename $rawGenome` -d ${PWD}/hic_data -D $juicer_path -z $rawGenome -y `basename $rawGenome`_${enzyme_type}.txt -p ${rawGenome}.sizes -s $enzyme_type -t 80
