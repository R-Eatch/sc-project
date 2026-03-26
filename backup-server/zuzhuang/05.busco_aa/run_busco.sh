
source /public/home/sunxuebo/software/miniconda3/etc/profile.d/conda.sh
export PATH=/public/home/sunxuebo/software/hmmer-3.3.2/bin:$PATH
export PATH=/public/home/fengchenguang/software/Augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/sunxuebo/software/busco-5.5.0/bin:$PATH

export AUGUSTUS_CONFIG_PATH=/public/home/sunxuebo/software/Augustus666/config





 
busco -i /data01/sunxuebo/project/zuzhuang/05.busco_aa/test.pep -o busco_out/ -l ../2.busco/mammalia_odb10 -m proteins --offline -f -c 30 --config /data01/sunxuebo/project/zuzhuang/1.3busco/config.ini




