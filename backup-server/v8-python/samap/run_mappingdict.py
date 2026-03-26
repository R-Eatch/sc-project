import re
import pandas as pd

# 输入pep文件和输出映射文件
pep_file = "R.pep"
mapping_file = "R.gene_mapping.csv"

# 存储映射关系
mapping = []

# 解析pep文件
with open(pep_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            # 从描述行中提取pep ID和基因ID
            match = re.search(r"^>(\S+).*gene:(\S+).*", line)
            if match:
                protein_id = match.group(1)  # 提取pep ID
                gene_id = match.group(2)     # 提取gene ID
                mapping.append((protein_id, gene_id))

# 转换为DataFrame并保存
df = pd.DataFrame(mapping, columns=["protein_id", "gene_id"])
df.to_csv(mapping_file, index=False)

print(f"Mapping extracted and saved to {mapping_file}")

