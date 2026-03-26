import subprocess
import pandas as pd
import os
from Bio import SeqIO

# 1. 使用Python方法从GFF文件中提取调控元件序列并保存为FASTA格式
def extract_sequences_from_gff(human_genome_fasta, gff_file, output_dir):
    genome = SeqIO.to_dict(SeqIO.parse(human_genome_fasta, "fasta"))
    gff_df = pd.read_csv(gff_file, sep='\t', header=None, comment='#', low_memory=False)
    gff_df.columns = ['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    
    output_fasta = f"{output_dir}/regulatory_elements.fa"
    with open(output_fasta, "w") as out_fasta:
        for _, row in gff_df.iterrows():
            chr_name = str(row['chr'])
            start = int(row['start'])  # 保持1-based坐标
            end = int(row['end'])
            
            if chr_name in genome:
                sequence = genome[chr_name].seq[start-1:end]  # 提取序列，注意GFF文件是1-based
                attributes = row['attributes']
                element_id = attributes.split(';')[0].split('=')[1]
                out_fasta.write(f">{element_id}\n{sequence}\n")
    
    return output_fasta

# 2. 使用BLAST将调控元件序列比对到蜜袋鼯基因组
def align_sequences_to_genome(query_file, sugar_glider_genome, output_dir, cores):
    db_path = f"{output_dir}/sugar_glider_genome"
    if not os.path.exists(f"{db_path}.nin"):
        subprocess.run(["makeblastdb", "-in", sugar_glider_genome, "-dbtype", "nucl", "-out", db_path])
    
    output_csv = f"{output_dir}/alignment.csv"
    subprocess.run([
        "blastn", "-query", query_file, "-db", db_path, 
        "-out", output_csv, "-outfmt", "6", "-num_threads", str(cores)
    ])
    return output_csv

# 3. 将BLAST结果转换为BED格式，并加入查询序列和分数信息
def convert_blast_to_bed(blast_output_csv, output_dir):
    # 读取BLAST输出CSV
    blast_df = pd.read_csv(blast_output_csv, header=None, delimiter='\t')  # 确保以制表符分隔
    
    # 检查列数是否正确，应该有12列
    if blast_df.shape[1] == 12:
        blast_df.columns = ['query', 'subject', 'identity', 'length', 'mismatch', 'gap_open', 
                            'q_start', 'q_end', 's_start', 's_end', 'e_value', 'bit_score']
    else:
        print(f"错误：BLAST输出的列数为{blast_df.shape[1]}，预期是12列。请检查输入文件。")
        return None
    
    # 创建BED文件格式数据，增加查询序列和分数
    bed_data = blast_df[['subject', 's_start', 's_end', 'query', 'bit_score']]
    bed_data['start'] = bed_data['s_start'] - 1  # 转换为0-based坐标
    bed_data['end'] = bed_data['s_end']
    
    bed_file = f"{output_dir}/alignment.bed"
    bed_data[['subject', 'start', 'end', 'query', 'bit_score']].to_csv(bed_file, sep='\t', header=False, index=False)
    
    return bed_file

# 4. 过滤BED文件：保留每个查询序列中分数最高的比对
def filter_bed_file(input_bed, output_bed):
    # 读取BED文件
    bed_df = pd.read_csv(input_bed, sep='\t', header=None)
    
    # 设置列名
    bed_df.columns = ['subject', 'start', 'end', 'query', 'bit_score']
    
    # 过滤：对于每个查询序列，保留分数最高的记录
    filtered_df = bed_df.loc[bed_df.groupby('query')['bit_score'].idxmax()]
    
    # 保存过滤后的BED文件
    filtered_df.to_csv(output_bed, sep='\t', header=False, index=False)
    
    print(f"过滤后的BED文件已保存到 {output_bed}")
    return output_bed

# 主函数
def main():
    human_genome_fasta = "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    gff_file = "homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20240230.gff"
    sugar_glider_genome = "/data01/sunxuebo/project/atac/data/sugarglider.fa"
    output_dir = "output"
    cores = 10  # 设置使用的核
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # 1. 提取调控元件序列
    regulatory_elements_fasta = extract_sequences_from_gff(human_genome_fasta, gff_file, output_dir)

    # 2. 将序列比对到蜜袋鼯基因组
    alignment_csv = align_sequences_to_genome(regulatory_elements_fasta, sugar_glider_genome, output_dir, cores)

    # 3. 将BLAST结果转换为BED格式
    bed_file = convert_blast_to_bed(alignment_csv, output_dir)

    # 4. 过滤BED文件：保留每个查询序列中分数最高的比对
    filtered_bed = filter_bed_file(bed_file, f"{output_dir}/filtered_alignment.bed")

    print(f"最终的过滤结果保存在 {filtered_bed}")

if __name__ == "__main__":
    main()

