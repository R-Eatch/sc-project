import pandas as pd

# 读取BED文件
bed_file = "output/filtered_alignment.bed"
bed_df = pd.read_csv(bed_file, sep='\t', header=None)

# 给BED文件列命名
bed_df.columns = ['chrom', 'start', 'end', 'name', 'score']

# 创建一个新的列'strand'来存储链方向
bed_df['strand'] = '.'

# 修正start和end的顺序（如果start大于end），并标记链方向
bed_df['strand'] = bed_df.apply(lambda row: '-' if row['start'] > row['end'] else '.', axis=1)

# 交换start和end的顺序，确保start小于end
bed_df['start'], bed_df['end'] = bed_df[['start', 'end']].min(axis=1), bed_df[['start', 'end']].max(axis=1)

# 输出修正后的BED文件
bed_df.to_csv("output/filtered_alignment_fixed.bed", sep='\t', header=False, index=False)

print("修正后的BED文件已保存为 output/filtered_alignment_fixed.bed")

