
input_file=./data/AG_TransLiG.fa
out_dir=./cds


/public/home/sunxuebo/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.LongOrfs -t ${input_file} -O $out_dir
/public/home/sunxuebo/software/TransDecoder-TransDecoder-v5.7.1/TransDecoder.Predict -t ${input_file} -O $out_dir

