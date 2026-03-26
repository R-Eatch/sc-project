#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/02.index.ref.sh.work/index_genome1
( time  /public/home/sunxuebo/software/NextPolish/bin/samtools faidx /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/input.genome.fasta )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/02.index.ref.sh.work/index_genome1/nextPolish.sh.done

