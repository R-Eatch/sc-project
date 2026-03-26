#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/03.map.ref.sh.work/map_genome2
( time  /public/home/sunxuebo/software/NextPolish/bin/minimap2 --split-prefix tmp -a -x asm20 -t 16 /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/input.genome.fasta /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/input.hifipart.001.fasta.gz|/public/home/sunxuebo/software/NextPolish/bin/samtools view --threads 5 -F 0x4 -b - |/public/home/sunxuebo/software/NextPolish/bin/samtools sort - -m 2g --threads 5 -o hifi.part001.sort.bam )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/03.map.ref.sh.work/map_genome2/nextPolish.sh.done

