samtools view -@ 32 -bS ../output_AG.sam > output_AG.bam
samtools view -@ 32 -bS ../output_MG.sam > output_MG.bam
samtools view -@ 32 -bS ../output_SG.sam > output_SG.bam


samtools merge -@ 32  HIC_combined.bam output_AG.bam output_MG.bam output_SG.bam



samtools sort -@ 32 HIC_combined.bam -o HIC_combined.sorted.bam

samtools index HIC_combined.sorted.bam

