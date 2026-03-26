#!/usr/bin/env python
# coding: utf-8

# In[3]:


import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
import scanpy as sc
import numpy as np
from venn import venn
n_topgenes=5000
re_find_topgenes = True


# In[1]:


# 定义文件路径和样本名称
file_paths = [
    "D:/111/R-MG_cleaned.h5ad",  # R-MG
    "D:/111/R-AG_cleaned.h5ad",  # R-AG
    "D:/111/S-MG_cleaned.h5ad",  # S-MG
    "D:/111/S-AG_cleaned.h5ad",  # S-AG
    "D:/111/M-MG_cleaned.h5ad",  # M-MG
    "D:/111/R-CG_cleaned.h5ad"   # R-CG
]
file_paths = [
    "../../R-MG/1.subset/R-MG_cleaned.h5ad",  # R-MG
    "../../R-AG/1.subset/R-AG_cleaned.h5ad",  # R-AG
    "../../S-MG/1.subset/S-MG_cleaned.h5ad",  # S-MG
    "../../S-AG/1.subset/S-AG_cleaned.h5ad",  # S-AG
    "../../M-MG/1.subset/M-MG_cleaned.h5ad",  # M-MG
    "../../R-CG/1.subset/R-CG_cleaned.h5ad"   # R-CG
]
sample_names = ["R-MG", "R-AG", "S-MG", "S-AG", "M-MG", "R-CG"]

# Step 1: 读取数据并提取 HVGs
def load_hvgs(file_paths, sample_names):
    hvgs_dict = {}
    for i, file_path in enumerate(file_paths):
        print(f"Loading sample: {sample_names[i]} from {file_path}...")
        adata = sc.read(file_path)
        if re_find_topgenes:
            adata.X=adata.layers['counts']
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.highly_variable_genes(adata, n_top_genes=n_topgenes)
            print(f're find hvgs: {n_topgenes}')
        hvgs = adata.var[adata.var['highly_variable']].index.tolist()  # 提取 HVGs
        hvgs_dict[sample_names[i]] = set(hvgs)  # 将 HVGs 存储为集合
        print(f"Sample {sample_names[i]}: {len(hvgs)} HVGs found.\n")
    return hvgs_dict

# Step 2: 计算交集数据
def calculate_intersections(hvgs_dict, sample_names):
    intersections = {}
    selected_combinations = []  # 存储感兴趣的组合

    # 两两交集
    for comb in combinations(sample_names, 2):
        selected_combinations.append(comb)

    # 三 MG 样本交集
    mg_samples = ["R-MG", "S-MG", "M-MG"]
    selected_combinations.append(tuple(mg_samples))

    # 所有六个样本的交集
    selected_combinations.append(tuple(sample_names))

    # 计算感兴趣的交集
    for comb in selected_combinations:
        intersect_set = set.intersection(*(hvgs_dict[sample] for sample in comb))
        intersections[comb] = len(intersect_set)

    return intersections

# Step 3: 绘制拼接图像
def plot_combined(intersections, sample_names, save_path="combined_plot.png"):
    # 准备数据
    combinations = list(intersections.keys())
    counts = list(intersections.values())
    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    combinations = [combinations[i] for i in sorted_indices]
    counts = [counts[i] for i in sorted_indices]

    # 创建图像布局
    fig = plt.figure(figsize=(20, 10))  # 增加宽度以拉长间距
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[2, 1], hspace=0.05)

    # 上方柱状图
    ax_bar = fig.add_subplot(grid[0, 0])
    bar_positions = range(len(combinations))
    bar_width = 0.6

    # 使用 viridis 配色
    cmap = plt.cm.get_cmap("viridis", len(combinations))
    colors = [cmap(i) for i in range(len(combinations))]

    bars = ax_bar.bar(bar_positions, counts, width=bar_width, color=colors)  # 使用颜色方案
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylabel("Number of HVGs")
    ax_bar.set_title("Intersection of HVGs Across Glandular Samples", fontsize=16)
    ax_bar.set_xlim(-0.5, len(combinations) - 0.5)  # 去掉左右空白

    # 在柱子上标注基因数量
    for bar, count in zip(bars, counts):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                    str(count), ha='center', fontsize=10)

    # 下方交集点线图
    ax_upset = fig.add_subplot(grid[1, 0])

    # 设置背景颜色交替
    for i in range(len(sample_names)):
        ax_upset.axhspan(i - 0.5, i + 0.5, color='lightgrey' if i % 2 == 0 else 'white', alpha=0.5)

    for i, comb in enumerate(combinations):
        for j, sample in enumerate(sample_names):
            if sample in comb:
                ax_upset.scatter(i, j, color='black', s=50)
        if len(comb) > 1:
            indices = [sample_names.index(sample) for sample in comb]
            ax_upset.plot([i, i], [min(indices), max(indices)], color='black', lw=1)

    ax_upset.set_xticks(bar_positions)  # 与柱状图共享 X 轴间距
    ax_upset.set_xticklabels([])
    ax_upset.set_yticks(range(len(sample_names)))
    ax_upset.set_yticklabels(sample_names)
    ax_upset.set_xlim(-0.5, len(combinations) - 0.5)  # 去掉左右空白
    ax_upset.set_ylim(-0.5, len(sample_names) - 0.5)

    # 保存图像
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.savefig("HVGs-for-all.pdf", dpi=400, bbox_inches='tight')
    plt.close()
    print(f"Combined plot saved as {save_path}")

# 主函数
def main(file_paths, sample_names):
    hvgs_dict = load_hvgs(file_paths, sample_names)
    intersections = calculate_intersections(hvgs_dict, sample_names)
    plot_combined(intersections, sample_names, save_path="HVGs-for-all.png")

# 运行主函数
main(file_paths, sample_names)


# In[ ]:


# 定义文件路径和样本名称
file_paths = [
    "../../R-MG/1.subset/R-MG_cleaned.h5ad",  # R-MG
    #"../../R-AG/1.subset/R-AG_cleaned.h5ad",  # R-AG
    "../../S-MG/1.subset/S-MG_cleaned.h5ad",  # S-MG
    #"../../S-AG/1.subset/S-AG_cleaned.h5ad",  # S-AG
    "../../M-MG/1.subset/M-MG_cleaned.h5ad",  # M-MG
   # "../../R-CG/1.subset/R-CG_cleaned.h5ad"   # R-CG
]
sample_names = ["R-MG", "S-MG", "M-MG"]

# Step 2: 计算任意两组的共有 HVGs 数量
def calculate_shared_hvgs(hvgs_dict, sample1, sample2):
    shared_hvgs = hvgs_dict[sample1] & hvgs_dict[sample2]  # 求交集
    print(f"Shared HVGs between {sample1} and {sample2}: {len(shared_hvgs)}")
    return shared_hvgs

# Step 3: 计算三个样本的共有 HVGs 数量
def calculate_all_shared_hvgs(hvgs_dict):
    # 计算三个集合的交集（三个样本的共有 HVGs）
    shared_hvgs_all = hvgs_dict["R-MG"] & hvgs_dict["S-MG"] & hvgs_dict["M-MG"]
    print(f"Shared HVGs across all samples: {len(shared_hvgs_all)}")
    return shared_hvgs_all

# Step 4: Venn 图可视化
def plot_venn(hvgs_dict, sample_names):
    # 使用 venn 绘制 3 集合的 Venn 图
    venn(hvgs_dict)
    plt.title("Venn Diagram of different species MG HVGs")
    plt.savefig("./3-MG-venn.png", dpi=300, bbox_inches='tight')
    plt.savefig("./3-MG-venn.pdf", dpi=300, bbox_inches='tight')
    plt.show()

# Step 5: 绘制柱状图展示不同交集的共有 HVGs 数量
def plot_bar_chart(shared_hvgs_dict):
    # shared_hvgs_dict 是一个字典，包含不同交集的 HVGs 数量
    categories = list(shared_hvgs_dict.keys())
    values = list(shared_hvgs_dict.values())

    # 绘制柱状图
    plt.figure(figsize=(8, 6))
    colors = ['#ff6666', '#6699cc', '#66cc66']  # 使用深色彩色柱状
    plt.bar(categories, values, color=colors, width=0.4)  # 设置柱宽
    plt.xlabel("Sample Groups", fontsize=12)
    plt.ylabel("Number of Shared HVGs", fontsize=12)
    plt.title("Shared HVGs Across Different Speices MG", fontsize=14)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig("./3-MG-barplot.png", dpi=300, bbox_inches='tight')
    plt.show()

# Step 6: 保存交集结果到 CSV
def save_shared_hvgs_to_csv(shared_hvgs_dict, hvgs_dict):
    all_shared_hvgs = {}
    for key, value in shared_hvgs_dict.items():
        if key == "R-MG & S-MG":
            all_shared_hvgs[key] = calculate_shared_hvgs(hvgs_dict, "R-MG", "S-MG")
        elif key == "R-MG & M-MG":
            all_shared_hvgs[key] = calculate_shared_hvgs(hvgs_dict, "R-MG", "M-MG")
        elif key == "S-MG & M-MG":
            all_shared_hvgs[key] = calculate_shared_hvgs(hvgs_dict, "S-MG", "M-MG")
        elif key == "All Samples (R-MG, S-MG, M-MG)":
            all_shared_hvgs[key] = calculate_all_shared_hvgs(hvgs_dict)

    # 保存到 CSV
    with open("shared_hvgs_MG.csv", "w") as f:
        f.write("Intersection,Shared HVGs\n")
        for group, hvgs in all_shared_hvgs.items():
            for hvg in hvgs:
                f.write(f"{group},{hvg}\n")
    print("Shared HVGs have been saved to shared_hvgs.csv")

# Step 7: 主函数
def main(file_paths, sample_names):
    # 加载 HVGs
    hvgs_dict = load_hvgs(file_paths, sample_names)

    # 打印并计算任意两组共有 HVGs 的个数
    for i in range(len(sample_names)):
        for j in range(i + 1, len(sample_names)):
            calculate_shared_hvgs(hvgs_dict, sample_names[i], sample_names[j])

    # 计算三个样本的共有 HVGs 数量
    shared_hvgs_all = calculate_all_shared_hvgs(hvgs_dict)

    # 打印共有 HVGs 数量
    print(f"\nShared HVGs across all samples: {len(shared_hvgs_all)}")

    # 创建包含共有 HVGs 数量的字典
    shared_hvgs_dict = {
        "R-MG & S-MG": len(calculate_shared_hvgs(hvgs_dict, "R-MG", "S-MG")),
        "R-MG & M-MG": len(calculate_shared_hvgs(hvgs_dict, "R-MG", "M-MG")),
        "S-MG & M-MG": len(calculate_shared_hvgs(hvgs_dict, "S-MG", "M-MG")),
        "All MGs": len(shared_hvgs_all)
    }

    # 绘制 Venn 图
    plot_venn(hvgs_dict, sample_names)

    # 绘制柱状图
    plot_bar_chart(shared_hvgs_dict)

    # 保存交集结果到 CSV
    save_shared_hvgs_to_csv(shared_hvgs_dict, hvgs_dict)

# 运行主函数
main(file_paths, sample_names)


# In[ ]:


file_paths = [
    "../../R-MG/1.subset/R-MG_cleaned.h5ad",  # R-MG
    "../../R-AG/1.subset/R-AG_cleaned.h5ad",  # R-AG
    "../../S-MG/1.subset/S-MG_cleaned.h5ad",  # S-MG
    "../../S-AG/1.subset/S-AG_cleaned.h5ad",  # S-AG
    "../../M-MG/1.subset/M-MG_cleaned.h5ad",  # M-MG
    "../../R-CG/1.subset/R-CG_cleaned.h5ad"   # R-CG
]
sample_names = ["R-MG", "R-AG", "S-MG", "S-AG", "M-MG", "R-CG"]


# Step 2: 计算交集并导出为 CSV
def calculate_and_export_intersections(hvgs_dict, sample_names, output_csv="pairwise_intersections.csv"):
    intersections = {}

    # 计算两两交集
    for sample1, sample2 in combinations(sample_names, 2):
        intersect_set = hvgs_dict[sample1] & hvgs_dict[sample2]
        intersection_name = f"{sample1} & {sample2}"
        intersections[intersection_name] = list(intersect_set)

    # 计算 M-MG, R-MG, S-MG 三交集
    three_intersection_set = hvgs_dict["M-MG"] & hvgs_dict["R-MG"] & hvgs_dict["S-MG"]
    intersections["M-MG & R-MG & S-MG"] = list(three_intersection_set)

    # 计算六样本交集
    six_intersection_set = set.intersection(*(hvgs_dict[sample] for sample in sample_names))
    intersections["All Samples Intersection"] = list(six_intersection_set)

    # 转换为 DataFrame 并导出 CSV
    max_length = max(len(genes) for genes in intersections.values())
    export_data = {name: genes + [""] * (max_length - len(genes)) for name, genes in intersections.items()}
    df = pd.DataFrame(export_data)
    df.to_csv(output_csv, index=False)
    print(f"Pairwise and additional intersections saved to {output_csv}")

# 主函数
def main(file_paths, sample_names):
    hvgs_dict = load_hvgs(file_paths, sample_names)
    calculate_and_export_intersections(hvgs_dict, sample_names, output_csv="HVGS-pairwise_intersections.csv")

# 运行主函数
main(file_paths, sample_names)


# In[ ]:


# 定义文件路径和样本名称
file_paths = [
    "../../R-MG/1.subset/R-MG_cleaned.h5ad",  # R-MG
    #"../../R-AG/1.subset/R-AG_cleaned.h5ad",  # R-AG
    "../../S-MG/1.subset/S-MG_cleaned.h5ad",  # S-MG
    #"../../S-AG/1.subset/S-AG_cleaned.h5ad",  # S-AG
    "../../M-MG/1.subset/M-MG_cleaned.h5ad",  # M-MG
    #"../../R-CG/1.subset/R-CG_cleaned.h5ad"   # R-CG
]
sample_names = ["R-MG", "S-MG", "M-MG"]

# Step 2: 计算各自独有的 HVGs
def calculate_unique_hvgs(hvgs_dict, sample_names):
    unique_hvgs_dict = {}
    for i, name in enumerate(sample_names):
        # 计算独有的 HVGs
        other_sets = [hvgs_dict[sample_names[j]] for j in range(len(sample_names)) if j != i]
        unique_hvgs = hvgs_dict[name] - set.union(*other_sets)
        unique_hvgs_dict[name] = unique_hvgs
        print(f"Unique HVGs for {name}: {len(unique_hvgs)}")
    return unique_hvgs_dict

# Step 3: 保存独有 HVGs 到 CSV
def save_unique_hvgs_to_csv(unique_hvgs_dict):
    all_data = []
    for sample, hvgs in unique_hvgs_dict.items():
        for hvg in hvgs:
            all_data.append([sample, hvg])
    
    df = pd.DataFrame(all_data, columns=["Sample", "HVG"])
    output_file = "unique_hvgs_MG.csv"
    df.to_csv(output_file, index=False)
    print(f"Unique HVGs have been saved to {output_file}")

# 主函数
def main(file_paths, sample_names):
    # 加载 HVGs
    hvgs_dict = load_hvgs(file_paths, sample_names)
    
    # 计算各自独有的 HVGs
    unique_hvgs_dict = calculate_unique_hvgs(hvgs_dict, sample_names)
    
    # 保存独有 HVGs 到 CSV
    save_unique_hvgs_to_csv(unique_hvgs_dict)

# 运行主函数
main(file_paths, sample_names)


# In[ ]:


# 文件路径设置
novel_gene_file = "/data02/sunxuebo/project/scrnaseq/no-mammal/mammalnewgene/mammalian_new_genes_final.csv"  # 包含新基因的 CSV
intersections_file = "HVGS-pairwise_intersections.csv"
output_csv = "HVGs_novel_gene_intersections.csv"

# Step 1: 读取新基因 CSV，并提取 Novel_Gene 列构成集合
novel_df = pd.read_csv(novel_gene_file)
# 过滤空值并转换为字符串，构成新基因集合
novel_genes = set(novel_df["mouse"].dropna().astype(str).tolist())
print(f"Total novel genes found: {len(novel_genes)}")

# Step 2: 读取交集 CSV
intersections_df = pd.read_csv(intersections_file)
print(f"Intersections CSV columns: {list(intersections_df.columns)}")

# Step 3: 针对每个交集类型，筛选出属于新基因的基因
novel_intersections = {}
for col in intersections_df.columns:
    # 取出该列数据，去除空值和空字符串
    genes = intersections_df[col].dropna().astype(str).tolist()
    genes = [g.strip() for g in genes if g.strip() != ""]
    # 筛选出新基因
    novel_inters = [g for g in genes if g in novel_genes]
    novel_intersections[col] = novel_inters
    print(f"Intersection '{col}': {len(novel_inters)} novel genes")

# Step 4: 找出所有交集类型中最长的列表长度，以便对齐格式
max_length = max(len(lst) for lst in novel_intersections.values())

# Step 5: 创建新 DataFrame，各列不足部分用空字符串填充
export_data = {}
for col, gene_list in novel_intersections.items():
    padded_list = gene_list + [""] * (max_length - len(gene_list))
    export_data[col] = padded_list

export_df = pd.DataFrame(export_data)

# Step 6: 导出 CSV
export_df.to_csv(output_csv, index=False)
print(f"Novel gene intersections saved to {output_csv}")


# In[ ]:




