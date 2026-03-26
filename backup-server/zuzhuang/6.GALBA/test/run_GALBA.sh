###env PATH!!!###


export PATH=/public/home/sunxuebo/software/GALBA-main/scripts:$PATH
export PATH=/public/home/sunxuebo/software/GALBA-main:$PATH
export PATH=/public/home/fengchenguang/software/Augustus/augustus-3.3.3/bin:$PATH
export AUGUSTUS_CONFIG_PATH=/public/home/sunxuebo/software/augustus-3.3.3/config
export AUGUSTUS_SCRIPTS_PATH=/public/home/sunxuebo/software/Augustus-master/scripts
export AUGUSTUS_BIN_PATH=/public/home/sunxuebo/software/augustus-3.3.3/bin
export PATH=/public/home/fengchenguang/software/Augustus/augustus-3.3.3/scripts:$PATH
export MINIPROT_PATH=/public/home/sunxuebo/software/miniprot-master  
export SCORER_PATH=/public/home/sunxuebo/software/miniprot-boundary-scorer-master
export MINIPROTHINT_PATH=/public/home/sunxuebo/software/miniprothint-master
export GENOMETHREADER_PATH=/public/home/sunxuebo/software/gth-1.7.3-Linux_x86_64-64bit/bin
export DIAMOND_PATH=/public/home/sunxuebo/software/diamond-linux
export TSEBRA_PATH=/public/home/sunxuebo/software/TSEBRA-main/bin

source /public/home/sunxuebo/software/miniconda3/etc/profile.d/conda.sh

conda activate GALBA_offline

##command test##
    galba.pl --species=midaiwu_contig --genome=contig_genome.fa \
       --prot_seq=Monodelphis_domestica_protein.fa \
      --threads 8
 
