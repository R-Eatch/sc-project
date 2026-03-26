#!/usr/bin/env python
"""
脚本功能：
1. 解析 GFF 文件，构建 protein_id 到 (gene, product) 的映射关系。
   GFF 文件中每行属性格式示例：
       gene=Prlr;product=prolactin receptor isoform X1;protein_id=XP_006520098.1
2. 读取小鼠蛋白 FASTA 文件，更新记录头为：
       >ProteinID GeneName Product [Mus musculus]
   其中 protein_id 为 FASTA 记录头（不含 “>”）；
   GeneName 和 Product 从 GFF 文件中获取。
3. 输出更新后的 FASTA 文件。
"""

import argparse
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="基于 GFF 文件更新小鼠蛋白 FASTA 文件的记录头格式"
    )
    parser.add_argument("-g", "--gff", required=True,
                        help="输入的 GFF 文件，包含属性字段：gene, product, protein_id")
    parser.add_argument("-f", "--fasta", required=True,
                        help="输入的小鼠蛋白 FASTA 文件")
    parser.add_argument("-o", "--output", default="updated_mouse_proteins.fasta",
                        help="输出更新后的 FASTA 文件 (default: updated_mouse_proteins.fasta)")
    return parser.parse_args()

def parse_gff(gff_file):
    """
    解析 GFF 文件，构建 mapping：
       protein_id -> (gene, product)
    适用于属性字段为：
       gene=Prlr;product=prolactin receptor isoform X1;protein_id=XP_006520098.1
    """
    mapping = {}
    with open(gff_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            # 如果 GFF 文件为标准格式，属性信息通常在最后一列，
            # 否则直接对整行进行解析
            parts = line.split("\t")
            if len(parts) >= 9:
                attributes = parts[8]
            else:
                attributes = line
            attr_dict = {}
            for item in attributes.split(";"):
                item = item.strip()
                if "=" in item:
                    key, value = item.split("=", 1)
                    attr_dict[key] = value
            # 判断是否包含所需的键
            if "protein_id" in attr_dict and "gene" in attr_dict and "product" in attr_dict:
                protein_id = attr_dict["protein_id"]
                gene = attr_dict["gene"]
                product = attr_dict["product"]
                mapping[protein_id] = (gene, product)
    return mapping

def update_fasta_headers(fasta_file, mapping, output_file):
    """
    遍历 FASTA 文件中的记录，根据 mapping 更新记录头为：
       >ProteinID GeneName Product [Mus musculus]
    若 protein_id 不在 mapping 中，则保持原记录头不变。
    """
    records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        protein_id = record.id  # 假设record.id为蛋白序列的ID，与GFF中的protein_id对应
        if protein_id in mapping:
            gene, product = mapping[protein_id]
            new_description = f"{protein_id} {gene} {product} [Mus musculus]"
            record.description = new_description
            # record.id 保持为 protein_id
        else:
            # 若未找到对应关系，可选择打印警告信息
            # print(f"Warning: {protein_id} 在GFF文件中未找到对应关系！")
            pass
        records.append(record)
    SeqIO.write(records, output_file, "fasta")
    print(f"更新后的 FASTA 文件已写入: {output_file}")

def main():
    args = parse_args()
    mapping = parse_gff(args.gff)
    print(f"从 GFF 文件中解析到 {len(mapping)} 条映射关系。")
    update_fasta_headers(args.fasta, mapping, args.output)

if __name__ == "__main__":
    main()

