import os
import re
import csv
import argparse
from Bio import SeqIO
from ete3 import Tree

# 哺乳动物物种定义
REPRESENTATIVE_MAMMALS = [
    "Homo sapien",
    "Mus musculus",
    "Bos taurus",
    "Macaca mulatta",
    "Ornithorhynchus anatinus",
    "Canis lupus familiaris",
    "Felis catus",
    "Monodelphis domestica",
    "Rattus norvegicus",
    "Equus caballus",
    "Oryctolagus cuniculus"
]

def parse_newick(file_path):
    """解析newick格式的树文件，去掉版本号并调整分支支持度位置"""
    with open(file_path, 'r') as file:
        newick_str = file.read().strip()
        
        # 去掉版本号（如 NP_001232435.1 -> NP_001232435）
        newick_str = re.sub(r'(\w+)\.\d+', r'\1', newick_str)
        
        # 确保支持度值放置在正确的位置
        newick_str = re.sub(r'(\d+\.\d+):(\d+\.\d+)', r':\2\1', newick_str)
    
    return newick_str

def load_sequences(fasta_files, mapping_file):
    """读取多个FASTA文件并返回ID和物种名的映射"""
    if os.path.exists(mapping_file):
        print(f"Loading species-to-sequence mapping from {mapping_file}...", flush=True)
        # 如果映射文件已存在，直接读取
        with open(mapping_file, 'r') as f:
            sequence_species = {line.split(',')[0]: line.strip().split(',')[1] for line in f.readlines()}
    else:
        print(f"Building species-to-sequence mapping...", flush=True)
        sequence_species = {}
        for fasta_file in fasta_files:
            for record in SeqIO.parse(fasta_file, "fasta"):
                # 物种名在描述中的[]中
                species_match = re.search(r'\[(.*?)\]', record.description)
                species = species_match.group(1) if species_match else "Unknown"
                # 去除版本号
                seq_id = re.sub(r'\.\d+$', '', record.id)
                sequence_species[seq_id] = species
        # 保存物种名与ID的映射文件
        with open(mapping_file, 'w') as f:
            for seq_id, species in sequence_species.items():
                f.write(f"{seq_id},{species}\n")
        print(f"Species-to-sequence mapping saved to {mapping_file}.", flush=True)

    return sequence_species

def classify_sequences(tree, sequence_species):
    """根据物种名称分类查询序列、哺乳动物序列和非哺乳动物序列"""
    query_sequences = set()
    mammalian_sequences = set()
    non_mammalian_sequences = set()

    for leaf in tree.iter_leaves():
        seq_id = leaf.name
        species = sequence_species.get(seq_id, "")
        
        if "Mus musculus" in species:
            query_sequences.add(seq_id)
        elif any(mammal in species for mammal in REPRESENTATIVE_MAMMALS):
            mammalian_sequences.add(seq_id)
        else:
            non_mammalian_sequences.add(seq_id)

    return query_sequences, mammalian_sequences, non_mammalian_sequences

def check_monophyly(query_sequences, mammalian_sequences, tree):
    """检查查询序列和哺乳动物序列是否单系"""
    try:
        common_ancestor = tree.get_common_ancestor(list(query_sequences) + list(mammalian_sequences))
        for leaf in common_ancestor.iter_leaves():
            if leaf.name not in mammalian_sequences and leaf.name not in query_sequences:
                return False
        return True
    except Exception as e:
        print(f"Error in check_monophyly: {e}", flush=True)
        return False

def process_tree(tree_file, sequence_species):
    """判断树是否为新基因Type1或Type2"""
    new_gene_type1 = []
    new_gene_type2 = []

    try:
        tree = Tree(parse_newick(tree_file), format=1)
    except Exception as e:
        print(f"Error parsing tree {tree_file}: {e}", flush=True)
        return [], []

    query_sequences_in_tree, mammalian_sequences_in_tree, non_mammalian_sequences_in_tree = classify_sequences(tree, sequence_species)

    if not non_mammalian_sequences_in_tree:
        # 如果树中所有序列都是哺乳动物序列和查询序列，认为是新基因Type1
        new_gene_type1.append(os.path.basename(tree_file).split('.')[0])  # 基因名（文件名）作为新基因Type1
    else:
        # 如果树中存在非哺乳动物序列，检查是否查询序列和哺乳动物序列单系
        if check_monophyly(query_sequences_in_tree, mammalian_sequences_in_tree, tree):
            new_gene_type2.append(os.path.basename(tree_file).split('.')[0])  # 基因名（文件名）作为新基因Type2

    return new_gene_type1, new_gene_type2

def save_to_csv(new_genes_type1, new_genes_type2, type1_file, type2_file):
    """保存新基因Type1和Type2的CSV文件"""
    with open(type1_file, 'w', newline='') as f1, open(type2_file, 'w', newline='') as f2:
        writer_type1 = csv.writer(f1)
        writer_type2 = csv.writer(f2)
        
        writer_type1.writerow(["GeneName"])
        writer_type2.writerow(["GeneName"])
        
        for gene in new_genes_type1:
            writer_type1.writerow([gene])
        
        for gene in new_genes_type2:
            writer_type2.writerow([gene])

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(description="分析系统发育树并判断新基因Type1和Type2")
    parser.add_argument("-q", "--query_fasta", required=True, help="查询序列FASTA文件")
    parser.add_argument("-d", "--database_fasta", required=True, help="数据库序列FASTA文件")
    parser.add_argument("-t", "--tree_dir", required=True, help="系统发育树文件目录")
    parser.add_argument("-o1", "--output_type1", default="new_genes_type1.csv", help="保存新基因Type1的CSV文件")
    parser.add_argument("-o2", "--output_type2", default="new_genes_type2.csv", help="保存新基因Type2的CSV文件")
    parser.add_argument("-m", "--mapping_file", default="species_to_sequence_mapping.csv", help="物种与序列ID映射文件")
    
    args = parser.parse_args()

    # 合并查询和数据库序列的FASTA文件，构建物种名映射
    sequence_species = load_sequences([args.query_fasta, args.database_fasta], args.mapping_file)
    
    # 用于存储新基因的列表
    new_genes_type1 = []
    new_genes_type2 = []

    # 遍历树文件并分析
    for filename in os.listdir(args.tree_dir):
        if filename.endswith(".nwk"):  # 假设树文件以.nwk结尾
            tree_file_path = os.path.join(args.tree_dir, filename)
            
            # 输出当前处理的树文件
            print(f"Processing {filename}...", flush=True)
            
            # 处理树文件，分类并判断新基因类型
            type1, type2 = process_tree(tree_file_path, sequence_species)
            new_genes_type1.extend(type1)
            new_genes_type2.extend(type2)

    # 保存新基因Type1和Type2的列表到CSV文件
    save_to_csv(new_genes_type1, new_genes_type2, args.output_type1, args.output_type2)
    print(f"Analysis complete. Results saved to {args.output_type1} and {args.output_type2}.", flush=True)

if __name__ == "__main__":
    main()

