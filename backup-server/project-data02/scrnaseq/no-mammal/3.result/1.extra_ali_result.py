#!/usr/bin/env python
"""
步骤 1：根据Diamond比对结果提取每个基因的比对结果并保存为FASTA文件。
1. 不再检查是否来源于哺乳动物。
2. 提取完成后，按照ID进行去重，确保每个基因的FASTA文件中ID唯一，并且每个FASTA文件仅包含该基因的比对结果。
"""
import os
from Bio import SeqIO
import argparse

def parse_diamond_results(diamond_file, query_fasta):
    """解析Diamond比对结果"""
    gene_hits = {}
    query_dict = {}
    
    # 读取查询FASTA文件，构建pep_id到基因名的映射
    with open(query_fasta, "r") as f:
        for record in SeqIO.parse(f, "fasta"):
            gene_name = record.description.split()[1]
            pep_id = record.id
            query_dict[pep_id] = gene_name

    # 解析Diamond比对结果
    with open(diamond_file, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 12:
                continue
            pep_id = parts[0]  # 第1列是pep ID
            sseqid = parts[1]  # 第2列是hit ID

            # 获取该pep_id对应的基因名
            gene_name = query_dict.get(pep_id)
            if gene_name:
                if gene_name not in gene_hits:
                    gene_hits[gene_name] = {}
                if pep_id not in gene_hits[gene_name]:
                    gene_hits[gene_name][pep_id] = set()  # 使用集合避免重复
                gene_hits[gene_name][pep_id].add(sseqid)

    return gene_hits

def extract_unique_sequences(gene_hits, query_fasta, database_fasta):
    """提取去重后的序列"""
    records = {}
    query_dict = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))
    db_dict = SeqIO.to_dict(SeqIO.parse(database_fasta, "fasta"))

    # 使用集合去重
    seen_ids = set()

    for gene_name, hits in gene_hits.items():
        gene_records = []  # 存储当前基因的FASTA记录

        # 提取查询序列
        for pep_id, hit_ids in hits.items():
            if pep_id in query_dict and pep_id not in seen_ids:
                gene_records.append(query_dict[pep_id])  # 添加查询序列
                seen_ids.add(pep_id)

        # 提取比对到的序列
        for pep_id, hit_ids in hits.items():
            for hit in hit_ids:
                if hit in db_dict and hit not in seen_ids:
                    gene_records.append(db_dict[hit])  # 添加比对到的序列
                    seen_ids.add(hit)

        # 只保存该基因的去重序列
        if gene_records:
            records[gene_name] = gene_records

    return records

def save_fasta_for_gene(gene_name, records):
    """保存FASTA文件"""
    output_dir = './data'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_file = os.path.join(output_dir, f"{gene_name}.fasta")
    with open(output_file, "w") as f:
        SeqIO.write(records, f, "fasta")

def main():
    parser = argparse.ArgumentParser(description="根据Diamond比对结果提取FASTA文件")
    parser.add_argument("-d", "--diamond", required=True, help="Diamond比对结果TSV文件")
    parser.add_argument("-q", "--query", required=True, help="查询pep文件（FASTA格式）")
    parser.add_argument("-b", "--database", required=True, help="数据库pep文件（FASTA格式）")
    args = parser.parse_args()

    # 解析Diamond结果，获得每个查询基因在比对上的命中集合
    gene_hits = parse_diamond_results(args.diamond, args.query)
    print(f"Diamond比对结果中共检测到 {len(gene_hits)} 个查询基因有命中。")

    # 提取去重后的序列并保存为基因名对应的FASTA文件
    gene_records = extract_unique_sequences(gene_hits, args.query, args.database)
    
    # 保存每个基因的去重FASTA文件
    for gene_name, records in gene_records.items():
        save_fasta_for_gene(gene_name, records)
    
    print(f"所有去重后的FASTA文件已保存到./data目录。")

if __name__ == "__main__":
    main()

