# OPTIONAL !!! if you run 3d-dna slowly, use split command to accelerate awk function
/public/home/zhuchenglong/software/seqkit/seqkit seq -w 80 ./sugar_glider.fa > ./input.fa


/public/home/zhengjiangmin/Software/3d-DNA/run-asm-pipeline.sh -r 2 ./input.fa hic_data/aligned/merged_nodups.txt
