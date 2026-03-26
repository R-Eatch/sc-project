#!/usr/bin/env python
"""
脚本功能：
1. 读取 CSV 文件中“genelist”列的基因名。
2. 解析小鼠蛋白 FASTA 文件（pep 文件），提取记录头中包含基因名的序列。
   假设 FASTA 记录头格式为：
      >ProteinID GeneName 描述信息 [Mus musculus]
3. 对于每个基因，提取所有对应的蛋白序列。
4. 将结果写入一个 FASTA 文件（默认文件名 query_proteins.fasta）。
"""

import csv
import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="从小鼠蛋白 FASTA 文件中提取 CSV 基因列表对应的蛋白序列"
    )
    parser.add_argument("-g", "--genelist", required=True,
                        help="CSV 文件，包含列 'genelist'（基因名称）")
    parser.add_argument("-f", "--fasta", required=True,
                        help="小鼠蛋白 FASTA 文件（pep 文件）")
    parser.add_argument("-o", "--output", default="query_proteins.fasta",
                        help="输出查询蛋白 FASTA 文件名 (default: query_proteins.fasta)")
    return parser.parse_args()

def read_gene_list(csv_file):
    """读取 CSV 中的基因名称，返回集合"""
    genes = set()
    with open(csv_file, 'r') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            gene = row['genelist'].strip()
            if gene:
                genes.add(gene)
    return genes

def extract_gene_from_header(header):
    """
    假设 header 格式为：
       >ProteinID GeneName 描述信息 [Mus musculus]
    则取第二个字段作为基因名称
    """
    parts = header.split()
    if len(parts) >= 2:
        return parts[1]
    else:
        return None

def select_all_proteins(records, target_genes):
    """
    遍历 FASTA 记录，对于目标基因的所有序列，
    提取每个基因对应的所有蛋白序列。
    返回字典：gene -> [SeqRecord]
    """
    gene_to_records = {}
    for record in records:
        gene = extract_gene_from_header(record.description)
        if gene in target_genes:
            if gene not in gene_to_records:
                gene_to_records[gene] = []
            gene_to_records[gene].append(record)
    return gene_to_records

def main():
    args = parse_args()
    target_genes = read_gene_list(args.genelist)
    print("加载目标基因数量：", len(target_genes))
    records = list(SeqIO.parse(args.fasta, "fasta"))
    print("FASTA 中总记录数：", len(records))
    gene_to_records = select_all_proteins(records, target_genes)
    print("找到对应蛋白序列的基因数量：", len(gene_to_records))
    # 写入所有找到的蛋白序列
    all_records = [record for records in gene_to_records.values() for record in records]
    SeqIO.write(all_records, args.output, "fasta")
    print("查询蛋白序列写入：", args.output)

if __name__ == "__main__":
    main()

