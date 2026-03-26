

#!/bin/bash

#/public/home/sunxuebo/software/yahs-main/juicer pre -a -o out_JBAT \
#    ../yahs.out.bin \
#    ../yahs.out_scaffolds_final.agp \
#    ../contig_genome.fa.fai
# -o out_JBAT表示输出文件名的前缀    

#JUICER=/public/home/sunxuebo/software/yahs-main/juicer_tools.2.20.00.jar
#asm_size=$(awk '{s+=$2} END{print s}' '../contig_genome.fa.fai')
#java -Xmx36G -jar $JUICER \
#    pre out_JBAT.txt out_JBAT.hic <(echo "assembly ${asm_size}")




/public/home/sunxuebo/software/yahs-main/juicer  pre -a -o out_JBAT ../yahs.out.bin ../yahs.out_scaffolds_final.agp ../contig_genome.fa.fai >out_JBAT.log 2>&1
(java -jar -Xmx32G /public/home/zhengjiangmin/Software/juicer-1.6/scripts/common/juicer_tools.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)

