#!/bin/bash
set -xveo pipefail
hostname
cd /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/01.db_split.sh.work/db_split1
( time  /public/home/sunxuebo/software/NextPolish/bin/seq_split -d /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir -m 500000000 -n 6 -i 0 -t 16 -f 1k -l 150k -s 355271384800 -p input.hifipart /data01/sunxuebo/project/zuzhuang/1.2nextpolish/hifi.fofn )
touch /data01/sunxuebo/project/zuzhuang/1.2nextpolish/polishdir/00.hifi_polish/01.db_split.sh.work/db_split1/nextPolish.sh.done

