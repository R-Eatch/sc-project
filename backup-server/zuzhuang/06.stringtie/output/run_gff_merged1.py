import re
import pandas as pd
from tqdm import tqdm

import pandas as pd

# 假设读取两个GFF文件
gff1_columns = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
gff1 = pd.read_csv(r"D:\111\sugar_glider_final.gff", sep='\t', names=gff1_columns, comment='#', header=None)
gff2 = pd.read_csv(r"D:\111\sugar_glider_final_exon.gff", sep='\t', names=gff1_columns, comment='#', header=None)
gff1['length'] = gff1['end'] - gff1['start']
df = gff2[gff2['attributes'].str.contains("reference_id", na=False)]
df = df[~df['type'].isin(['transcript'])]
df['length'] = df['end'] - df['start']
df['ID'] = df['attributes'].str.extract(r'reference_id\s*\"([^\"]+)\"')
df0 = gff1
df0['ID'] = gff1['attributes'].str.extract(r'ID=([^\;]+)')
merged_gff = pd.concat([df0, df], ignore_index=True)
merged_gff.sort_values(by=['ID', 'start'], inplace=True)
from tqdm import tqdm  # 用于进度条

# 假设 df 是你的 DataFrame，包含 'ID', 'type', 'length' 列

# 用 groupby 减少切片操作
grouped = merged_gff.groupby('ID')

# 初始化一个空列表，用于存储过滤后的 DataFrame
filtered_list = []

# 使用 tqdm 包裹循环，显示进度条
for id_value, id_df in tqdm(grouped, desc="Processing IDs"):
    # 提取 type 为 mRNA 的行并获取其 length
    mrna_length = id_df.loc[id_df['type'] == 'mRNA', 'length'].values

    if len(mrna_length) > 0:
        mrna_length = mrna_length[0]  # 假设每个 ID 只有一个 mRNA

        # 使用条件过滤，保留 type 为 exon 且 length 不等于 mRNA length 的行
        exon_filtered = id_df[~((id_df['type'] == 'exon') & (id_df['length'] == mrna_length))]

        # 将过滤后的结果添加到列表中
        filtered_list.append(exon_filtered)
    else:
        # 如果没有 mRNA 行，直接添加当前组的所有行
        filtered_list.append(id_df)

# 在循环结束后，将所有结果一次性拼接为 DataFrame
filtered_df = pd.concat(filtered_list, ignore_index=True)
import re

# 用 groupby 按照 ID 分组
grouped = filtered_df.groupby('ID')

# 定义一个空列表存储处理后的结果
updated_list = []

# 使用 tqdm 显示进度条
for id_value, id_df in tqdm(grouped, desc="Processing IDs"):
    # 提取 type=gene 的行，获取 attributes 列中的 Name=genename 信息
    gene_row = id_df[id_df['type'] == 'gene']

    if not gene_row.empty:
        # 使用正则表达式提取 Name=genename 部分的值
        gene_name_match = re.search(r'Name=([^;]+);', gene_row['attributes'].values[0])
        if gene_name_match:
            gene_name = gene_name_match.group(1)  # 提取到的基因名，如 "ACTA2"
            # 构建 gene_id "ACTA2" 字符串
            new_gene_id = f'gene_id "{gene_name}"'

            # 遍历该 ID 下所有 type=exon 的行，并替换 gene_id 信息
            id_df['attributes'] = id_df.apply(
                lambda row: re.sub(r'gene_id\s*"[^"]+"', new_gene_id, row['attributes'])
                if row['type'] == 'exon' else row['attributes'], axis=1
            )

    # 将处理后的数据添加到列表中
    updated_list.append(id_df)

# 最后，将所有的处理结果拼接成一个 DataFrame
updated_df = pd.concat(updated_list, ignore_index=True)

# 指定要保存的列
columns_to_save = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

# 过滤出只包含指定列的 DataFrame
final_df = updated_df[columns_to_save]

# 将 DataFrame 保存为 GFF 文件，不保存索引，并且不包含 header
final_df.to_csv("sugar_glider_final_have_exon.gff", sep='\t', index=False, header=False)

