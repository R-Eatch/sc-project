#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/05.polish.ref.sh.work/polish_genome1
( time  /public/home/sunxuebo/software/miniconda3/bin/python /public/home/sunxuebo/software/NextPolish/lib/nextpolish2.py -sp -p 16 -g /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/input.genome.fasta -b /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/input.genome.fasta.blc -i 0 -l /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/hifi.sort.bam.list -r hifi -o genome.nextpolish.part000.fasta )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/05.polish.ref.sh.work/polish_genome1/nextPolish.sh.done

