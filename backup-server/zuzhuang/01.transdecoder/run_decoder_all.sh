#!/bin/bash

data_dir="./data"
out_dir="./cds"
transdecoder_path="/public/home/sunxuebo/software/TransDecoder-TransDecoder-v5.7.1"

for input_file in $data_dir/*TransLiG.fa; do
    organ_name=$(basename "$input_file" | cut -d'_' -f1)
    pbs_script="decode_${organ_name}.pbs"

    echo "#!/bin/bash" > $pbs_script
    echo "#PBS -N decode_${organ_name}" >> $pbs_script
    echo "#PBS -o ./" >> $pbs_script
    echo "#PBS -e ./" >> $pbs_script
    echo "#PBS -q more" >> $pbs_script
    echo "#PBS -l nodes=1:ppn=2" >> $pbs_script
    echo "#PBS -r y" >> $pbs_script
    echo "" >> $pbs_script
    echo "cd $PWD" >> $pbs_script
    echo "$transdecoder_path/TransDecoder.LongOrfs -t ${input_file} -O $out_dir" >> $pbs_script
    echo "$transdecoder_path/TransDecoder.Predict -t ${input_file} -O $out_dir" >> $pbs_script

    qsub $pbs_script
done

