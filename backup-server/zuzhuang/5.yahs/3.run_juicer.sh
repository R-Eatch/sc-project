

#(/public/home/sunxuebo/software/yahs-main/juicer pre yahs.out.bin yahs.out_scaffolds_final.agp  contig_genome.fa.fai | sort -k2,2d -k6,6d -T ./ --parallel=10 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

#/public/home/sunxuebo/software/bioawk-1.0/bioawk -c fastx '{print $name, length($seq)}' yahs.out_scaffolds_final.fa > yahs.out_scaffolds_final.sizes




(java -jar -Xmx200G /public/home/sunxuebo/software/yahs-main/juicer_tools.2.20.00.jar  pre alignments_sorted.txt out.hic.part yahs.out_scaffolds_final.sizes) && (mv out.hic.part out.hic)
