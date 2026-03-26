import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

# 定义哺乳动物的物种名列表
mammal_species = [
    "Homo sapien",
    "Bos taurus",
    "Macaca mulatta",
    "Ornithorhynchus anatinus",
    "Canis lupus familiaris",
    "Felis catus",
    "Monodelphis domestica",
    "Rattus norvegicus",
    "Equus caballus",
    "Oryctolagus cuniculus"
    # 在这里添加更多的哺乳动物物种名
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

def analyze_genes(csv_file, data_dir, mammal_species, result_dir):
    """读取基因列表并分析每个基因的FASTA文件"""
    # 读取CSV文件，假设CSV文件中有基因名称列
    gene_data = pd.read_csv(csv_file)
    
    # 用于存储结果的列表
    results = []

    # 用于存储每个哺乳动物种类数下的基因名称
    species_genes_map = defaultdict(list)

    # 遍历基因列表
    for gene in gene_data['GeneName']:
        fasta_file = os.path.join(data_dir, f"{gene}.fasta")
        if not os.path.exists(fasta_file):
            print(f"Warning: {fasta_file} not found. Skipping gene {gene}.")
            continue

        # 提取物种信息
        species_list, species_set = extract_species_from_fasta(fasta_file)

        # 计算基因涉及的哺乳动物物种数量
        mammal_species_in_gene = [species for species in species_list if species in mammal_species]
        num_sequences = len(species_list)
        num_mammals = len(mammal_species_in_gene)
        num_mammal_species = len(set(mammal_species_in_gene))  # 唯一的哺乳动物物种数量

        # 存储结果
        results.append([gene, num_sequences, num_mammals, num_mammal_species, len(species_set)])

        # 将基因按哺乳动物种类数进行分类
        species_genes_map[num_mammal_species].append(gene)

    # 检查结果目录是否存在，不存在则创建
    os.makedirs(result_dir, exist_ok=True)

    # 将结果存入CSV文件
    results_df = pd.DataFrame(results, columns=['Gene', 'Total Sequences', 'Mammal Sequences', 'Mammal Species Count', 'Total Species Count'])
    results_df.to_csv(os.path.join(result_dir, 'gene_analysis_results.csv'), index=False)

    # 导出每个哺乳动物种类数下的基因名称
    species_genes_df = pd.DataFrame({k: pd.Series(v) for k, v in species_genes_map.items()})
    species_genes_df.to_csv(os.path.join(result_dir, 'genes_by_mammal_species_count.csv'), index=False)

    return results_df

def plot_histogram(data, result_dir):
    """绘制直方图，展示不同哺乳动物物种数量的基因分布"""
    # 统计每个基因的哺乳动物物种数量
    species_count = data['Mammal Species Count'].value_counts().sort_index()

    # 绘制直方图
    plt.figure(figsize=(8, 6))
    ax = species_count.plot(kind='bar', color='skyblue', edgecolor='black')

    plt.title('Distribution of Genes by Number of Mammal Species Involved')
    plt.xlabel('Number of Mammal Species')
    plt.ylabel('Gene Count')
    plt.xticks(rotation=0)

    # 在每个柱子顶部显示基因个数
    for i, v in enumerate(species_count):
        ax.text(i, v + 0.1, str(v), ha='center', va='bottom', fontweight='bold')

    # 保存图像到指定目录
    plt.savefig(os.path.join(result_dir, 'mammal_species_histogram.png'))
    plt.show()

# 主程序
if __name__ == "__main__":
    # 设置要读取的基因CSV文件路径和数据目录
    csv_file = 'new_genes_type1.csv'  # 替换为实际的基因列表文件路径
    data_dir = './data/'  # 替换为FASTA文件所在目录路径
    result_dir = './4.result/'  # 设置结果导出目录

    # 进行基因分析
    gene_analysis_results = analyze_genes(csv_file, data_dir, mammal_species, result_dir)

    # 绘制并保存直方图
    plot_histogram(gene_analysis_results, result_dir)

