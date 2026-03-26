#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/01.merge.ref.sh.work/merge_ref1
( time  cat /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/04.polish.ref.sh.work/polish_genome1/genome.nextpolish.part000.fasta /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/04.polish.ref.sh.work/polish_genome2/genome.nextpolish.part001.fasta /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/04.polish.ref.sh.work/polish_genome3/genome.nextpolish.part002.fasta > /data01/sunxuebo/project/zuzhuang/1.2nextpolish/./polishdir/01.hifi_polish/input.genome.fasta )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/01.hifi_polish/01.merge.ref.sh.work/merge_ref1/merge_ref.sh.done

