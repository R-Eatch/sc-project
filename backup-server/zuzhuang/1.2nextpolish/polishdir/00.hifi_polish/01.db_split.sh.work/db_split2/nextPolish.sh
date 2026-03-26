#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/01.db_split.sh.work/db_split2
( time  /public/home/sunxuebo/software/NextPolish/bin/samtools faidx /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/input.genome.fasta )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/01.db_split.sh.work/db_split2/nextPolish.sh.done

