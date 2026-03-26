#!/bin/bash
# 1. 解压缩当前目录下所有的 zip 文件
echo "解压缩当前目录下所有的 zip 文件..."
for zip_file in *.zip; do
  if [ -f "$zip_file" ]; then
    echo "正在解压：$zip_file"
    unzip -o "$zip_file"
  fi
done

# 2. 查找 protein.faa 文件并合并
echo "查找当前目录及子目录中的 protein.faa 文件..."
faa_files=$(find . -type f -name "protein.faa")
num_files=$(echo "$faa_files" | wc -l)
echo "找到 $num_files 个 protein.faa 文件。"
if [ "$num_files" -ne 20 ]; then
  echo "警告：预期应为 20 个 protein.faa 文件，但实际找到 $num_files 个。"
fi

# 3. 合并所有 protein.faa 文件为一个 FASTA 文件
output_file="merged_fasta_pep_no_mammals.fasta"
echo "合并所有 protein.faa 文件到 $output_file ..."
find . -type f -name "protein.faa" -exec cat {} + > "$output_file"
echo "合并完成。输出文件为：$output_file"

