import pandas as pd

# 输入文件路径
blast_file = "ms_to_rb.txt"
output_file = "ms_to_rb1.txt"
query_mapping_file = "../../M.gene_mapping.csv"
subject_mapping_file = "../../R.gene_mapping.csv"

# 是否替换query_protein和subject_protein
replace_query = True  # 设置为 False 禁用 query_protein 的替换
replace_subject = True  # 设置为 False 禁用 subject_protein 的替换

# 读取BLAST结果
blast_results = pd.read_csv(blast_file, sep="\t", header=None)
blast_results.columns = ["query_protein", "subject_protein", "percent_identity", "alignment_length",
                         "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]

# 读取映射文件并设置索引以加快查询
if replace_query:
    query_mapping = pd.read_csv(query_mapping_file)
    query_mapping = query_mapping.set_index("protein_id")  # 设置索引为 protein_id
else:
    query_mapping = None

if replace_subject:
    subject_mapping = pd.read_csv(subject_mapping_file)
    subject_mapping = subject_mapping.set_index("protein_id")  # 设置索引为 protein_id
else:
    subject_mapping = None

# 替换 query_protein 列
if replace_query:
    blast_results["query_gene"] = blast_results["query_protein"].map(query_mapping["gene_id"])
else:
    blast_results["query_gene"] = blast_results["query_protein"]  # 保留原始列

# 替换 subject_protein 列
if replace_subject:
    blast_results["subject_gene"] = blast_results["subject_protein"].map(subject_mapping["gene_id"])
else:
    blast_results["subject_gene"] = blast_results["subject_protein"]  # 保留原始列

# 删除原始 protein 列，保留 gene 列
final_results = blast_results[["query_gene", "subject_gene", "percent_identity", "alignment_length",
                                "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]]

# 保存到新文件
final_results.to_csv(output_file, sep="\t", index=False, header=False)

print(f"Updated BLAST results saved to {output_file}")

