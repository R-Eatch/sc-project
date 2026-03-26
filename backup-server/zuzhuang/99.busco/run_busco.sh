
source /public/home/sunxuebo/software/miniconda3/etc/profile.d/conda.sh
export PATH=/public/home/sunxuebo/software/hmmer-3.3.2/bin:$PATH
export PATH=/public/home/fengchenguang/software/Augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/sunxuebo/software/busco-5.5.0/bin:$PATH

cd /data01/sunxuebo/project/zuzhuang/99.busco
export AUGUSTUS_CONFIG_PATH=/public/home/sunxuebo/software/Augustus666/config






 busco -i yahs.out_scaffolds_final.fa  -o busco_out/ -l ../2.busco/mammalia_odb10 -m genome --offline -r --augustus -c 30 --config /data01/sunxuebo/project/zuzhuang/1.3busco/config.ini
 




