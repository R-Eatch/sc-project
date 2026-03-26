#!/bin/bash

# 目标目录路径
target_directory="/data01/sunxuebo/project/zuzhuang/01.transdecoder/data"

# 确保目标目录存在
mkdir -p "$target_directory"

# 循环遍历每个文件路径
for path in \
./AG/TransLiG_Out_Dir/TransLiG.fa \
./brain/TransLiG_Out_Dir/TransLiG.fa \
./heart/TransLiG_Out_Dir/TransLiG.fa \
./intestines/TransLiG_Out_Dir/TransLiG.fa \
./kidney/TransLiG_Out_Dir/TransLiG.fa \
./liver-1/TransLiG_Out_Dir/TransLiG.fa \
./lung/TransLiG_Out_Dir/TransLiG.fa \
./marrow/TransLiG_Out_Dir/TransLiG.fa \
./MG/TransLiG_Out_Dir/TransLiG.fa \
./muscle/TransLiG_Out_Dir/TransLiG.fa \
./ovary/TransLiG_Out_Dir/TransLiG.fa \
./SG/TransLiG_Out_Dir/TransLiG.fa \
./skin/TransLiG_Out_Dir/TransLiG.fa \
./spleen/TransLiG_Out_Dir/TransLiG.fa \
./stomach/TransLiG_Out_Dir/TransLiG.fa \
./uterus/TransLiG_Out_Dir/TransLiG.fa \
./pancrease/TransLiG_Out_Dir/TransLiG.fa
do
  # 从文件路径中提取器官名
  organ_name=$(basename $(dirname $(dirname $path)))
  
  # 构建新的文件名
  new_file_name="${organ_name}_TransLiG.fa"
  
  # 复制并重命名文件到目标目录
  cp "$path" "$target_directory/$new_file_name"
done

echo "所有文件已复制到 $target_directory"

