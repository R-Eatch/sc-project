#!/bin/bash

data_dir="../proteins"
out_dir="./"

for input_file in $data_dir/*proteins.fasta; do
    organ_name=$(basename "$input_file" | cut -d'_' -f1)
    pbs_script="predict_${organ_name}.pbs"

    echo "#!/bin/bash" > $pbs_script
    echo "#PBS -N predict_${organ_name}" >> $pbs_script
    echo "#PBS -o ./" >> $pbs_script
    echo "#PBS -e ./" >> $pbs_script
    echo "#PBS -q more" >> $pbs_script
    echo "#PBS -l nodes=1:ppn=20" >> $pbs_script
    echo "#PBS -r y" >> $pbs_script
    echo "" >> $pbs_script
    echo "cd $PWD" >> $pbs_script
    echo "source /public/software/profile.d/compiler_gnu-7.2.0.sh
          export CC=/public/software/compiler/gnu/7.2.0/bin/gcc
          export CXX=/public/software/compiler/gnu/7.2.0/bin/g++"  >> $pbs_script
    echo "/data01/zhuchenglong/software/python3/Python-3.9.16_build/bin/python3 /data01/wangkun/projects/gene_prediction/version4/test/miniprotGeneAnnotation.py -g sugar_glider.fa -p ${data_dir}/${organ_name}_proteins.fasta -t 20 -o ./${organ_name}.gff" >> $pbs_script

done

