import os
import pandas as pd
from collections import defaultdict

# 定义需要同时包含的五个物种的名称
required_species = [
    "Mus musculus",  # 小鼠
    "Oryctolagus cuniculus",  # 兔
    "Monodelphis domestica",  # 短尾负鼠
    "Ornithorhynchus anatinus",  # 鸭嘴兽
    "Homo sapiens"  # 人类
]

def extract_species_from_fasta(fasta_file):
    """从FASTA文件中提取物种名称"""
    species_list = []
    species_set = set()  # 用于统计物种种类数
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):  # 只处理以 '>' 开头的抬头行
                # 提取物种名称，格式为 [Homo sapiens] 等
                species_name = line.split('[')[-1].split(']')[0].strip()
                species_list.append(species_name)
                species_set.add(species_name)
    return species_list, species_set

def analyze_genes(csv_file, data_dir, required_species, result_dir):
    """读取基因列表并筛选包含特定五个物种的基因"""
    # 读取CSV文件，假设CSV文件中有基因名称列
    gene_data = pd.read_csv(csv_file)
    
    # 用于存储符合条件的结果的列表
    results = []

    # 遍历基因列表
    for gene in gene_data['GeneName']:
        fasta_file = os.path.join(data_dir, f"{gene}.fasta")
        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} not found. Skipping gene {gene}.")
            continue

        # 提取物种信息
        species_list, species_set = extract_species_from_fasta(fasta_file)

        # 计算该基因是否至少包含指定的五个物种
        if all(species in species_set for species in required_species):
            num_sequences = len(species_list)
            num_mammals = len([species for species in species_list if species in required_species])
            num_mammal_species = len(set(required_species).intersection(species_set))  # 只统计符合五个物种要求的

            # 存储结果
            results.append([gene, num_sequences, num_mammals, num_mammal_species, len(species_set)])

    # 检查结果目录是否存在，不存在则创建
    os.makedirs(result_dir, exist_ok=True)

    # 将符合条件的结果存入CSV文件
    results_df = pd.DataFrame(results, columns=['Gene', 'Total Sequences', 'Mammal Sequences', 'Mammal Species Count', 'Total Species Count'])
    results_df.to_csv(os.path.join(result_dir, 'filtered_genes_results.csv'), index=False)

    return results_df

# 主程序
if __name__ == "__main__":
    # 设置要读取的基因CSV文件路径和数据目录
    csv_file = 'new_genes_type1.csv'  # 替换为实际的基因列表文件路径
    data_dir = './data/'  # 替换为FASTA文件所在目录路径
    result_dir = './4.result/'  # 设置结果导出目录

    # 进行基因分析
    filtered_genes_results = analyze_genes(csv_file, data_dir, required_species, result_dir)

    # 显示结果（可选）
    print(filtered_genes_results)

