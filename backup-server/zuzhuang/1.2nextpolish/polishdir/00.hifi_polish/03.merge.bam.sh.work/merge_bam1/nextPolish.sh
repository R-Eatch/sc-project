#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/03.merge.bam.sh.work/merge_bam1
( time  /public/home/sunxuebo/software/NextPolish/bin/samtools merge -f -b /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/00.hifi_polish/hifi.sort.bam.list --threads 8 /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/00.hifi_polish/hifi.sort.bam )
( time  /public/home/sunxuebo/software/NextPolish/bin/samtools index -@ 16 /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/00.hifi_polish/hifi.sort.bam )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/03.merge.bam.sh.work/merge_bam1/nextPolish.sh.done

