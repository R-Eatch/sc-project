
source /public/home/sunxuebo/software/miniconda3/etc/profile.d/conda.sh
export PATH=/public/home/sunxuebo/software/hmmer-3.3.2/bin:$PATH
export PATH=/public/home/fengchenguang/software/Augustus/augustus-3.3.3/bin:$PATH
export PATH=/public/home/sunxuebo/software/busco-5.5.0/bin:$PATH

cd /data01/sunxuebo/project/zuzhuang/2.busco
export AUGUSTUS_CONFIG_PATH=/public/home/sunxuebo/software/Augustus666/config






 busco -i ./p_ctg.fa -o busco_out/p_ctg -l ./mammalia_odb10 -m genome --offline -f --augustus -c 50 --config /data01/sunxuebo/project/zuzhuang/2.busco/config.ini
 
# busco -i ./hap1_p_ctg.fa -o busco_out/hap1_p_ctg -l ./mammalia_odb10 -m genome --offline -f  --augustus
# busco -i ./hap2_p_ctg.fa -o busco_out/hap2_p_ctg -l ./mammalia_odb10 -m genome --offline -f  --augustus



