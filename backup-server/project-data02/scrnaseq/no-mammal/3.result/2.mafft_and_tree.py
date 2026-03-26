#!/usr/bin/env python
"""
步骤 2：对当前目录 `data` 文件夹中的FASTA文件进行MAFFT多序列比对，并使用FastTree构建系统发育树。
只保留树文件，命名格式为基因名。
"""
import os
import subprocess
from Bio import SeqIO
import argparse

def run_mafft(input_fasta, output_aln, threads="8"):
    """调用MAFFT进行多序列比对，支持多线程"""
    cmd = ["mafft", "--thread", threads, input_fasta]
    with open(output_aln, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

def run_fasttree(input_aln, output_tree):
    """调用FastTree构建系统发育树"""
    cmd = ["fasttree", input_aln]
    with open(output_tree, "w") as out_f:
        subprocess.run(cmd, stdout=out_f, check=True)

def main():
    parser = argparse.ArgumentParser(description="基于当前目录 `data` 文件夹中的FASTA文件进行MAFFT比对和FastTree树构建")
    parser.add_argument("-d", "--data_dir", default="data", help="包含FASTA文件的目录")
    parser.add_argument("-o", "--output_dir", default="output_trees", help="保存树文件的目录")
    parser.add_argument("--threads", default="8", help="MAFFT使用的线程数")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # 遍历data文件夹中的所有FASTA文件
    for filename in os.listdir(args.data_dir):
        if filename.endswith(".fasta"):
            gene_name = filename.split(".")[0]  # 假设文件名就是基因名
            input_fasta = os.path.join(args.data_dir, filename)
            aln_file = os.path.join(args.output_dir, f"{gene_name}.aln")
            tree_file = os.path.join(args.output_dir, f"{gene_name}.nwk")

            try:
                # 进行MAFFT比对
                run_mafft(input_fasta, aln_file, threads=args.threads)
                # 使用FastTree构建系统发育树
                run_fasttree(aln_file, tree_file)
                print(f"基因 {gene_name} 的树文件已保存到 {tree_file}")
                # 删除对齐文件，只保留树文件
                os.remove(aln_file)
            except Exception as e:
                print(f"处理基因 {gene_name} 时出错: {e}")

if __name__ == "__main__":
    main()

