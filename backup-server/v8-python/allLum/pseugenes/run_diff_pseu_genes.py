#!/usr/bin/env python
# coding: utf-8

# In[12]:


from __future__ import annotations
import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Sequence

# 定义文件路径和样本名称（注意CSV文件的命名格式）
file_paths = [
        "/data01/sunxuebo/project/scrnaseq/v8-python/R-MG/9.stage_DEG/R-MG_HR-SEC_pseudotime_genes.csv",
        "/data01/sunxuebo/project/scrnaseq/v8-python/R-AG/9.stage_DEG/R-AG_HR-SEC_pseudotime_genes.csv",
        "/data01/sunxuebo/project/scrnaseq/v8-python/S-MG/9.stage_DEG/S-MG_HR-SEC_pseudotime_genes.csv",
        "/data01/sunxuebo/project/scrnaseq/v8-python/S-AG/9.stage_DEG/S-AG_HR-SEC_pseudotime_genes.csv",
        "/data01/sunxuebo/project/scrnaseq/v8-python/M-MG/9.stage_DEG/M-MG_HR-SEC_pseudotime_genes.csv",
        "/data01/sunxuebo/project/scrnaseq/v8-python/R-CG/9.stage_DEG/R-CG_HR-SEC_pseudotime_genes.csv",
]
sample_names = ["R-MG", "R-AG", "S-MG", "S-AG", "M-MG", "R-CG"]
# ---------------- configurable parameters ----------------
TOP_N: int = 300          # number of top genes kept *per stage*
PVAL_CUTOFF: float = 0.01
LFC_CUTOFF: float = 0.25
# ---------------------------------------------------------
# Step 1: 读取数据并提取基因列表（取CSV文件第一列）
# def load_genelists(file_paths, sample_names):
#     gene_dict = {}
#     for i, file_path in enumerate(file_paths):
#         print(f"Loading sample: {sample_names[i]} from {file_path}...")
#         df = pd.read_csv(file_path)
#         # 假设CSV文件的第一列为基因名称
#         gene_list = df.iloc[:, 0].tolist()
#         gene_dict[sample_names[i]] = set(gene_list)
#         print(f"Sample {sample_names[i]}: {len(gene_list)} genes found.\n")
#     return gene_dict
def _select_top_by_stage(df: pd.DataFrame,
                         stage_col: str,
                         gene_col: str,
                         pval_cutoff: float,
                         lfc_cutoff: float,
                         top_n: int) -> Sequence[str]:
    """Return union of top_n genes per stage in *df* (filtered by sig)."""
    stage_top: Dict[str, List[str]] = {}
    for st, sub in df.groupby(stage_col):
        sig = sub[(sub["pvals_adj"] < pval_cutoff) & (sub["scores"] > lfc_cutoff)]
        top_genes = (
            sig.sort_values("scores", ascending=False)[gene_col]
            .head(top_n)
            .tolist()
        )
        stage_top[st] = top_genes
    return pd.unique(np.concatenate(list(stage_top.values())))


def load_genelists(file_paths: List[str],
                   sample_names: List[str],
                   top_n: int = TOP_N,
                   pval_cutoff: float = PVAL_CUTOFF,
                   lfc_cutoff: float = LFC_CUTOFF) -> Dict[str, Set[str]]:
    """Read the six CSVs and return {sample_name: set(genes)}.

    *file_paths* and *sample_names* **must** align by index.
    """
    if len(file_paths) != len(sample_names):
        raise ValueError("file_paths and sample_names must have the same length")

    gene_dict: Dict[str, Set[str]] = {}
    for path, sample in zip(file_paths, sample_names):
        csv_path = Path(path)
        if not csv_path.exists():
            raise FileNotFoundError(csv_path)

        print(f"Processing {sample} …")
        df = pd.read_csv(csv_path)

        # Detect column names
        stage_col = df.columns[0]                 # first column holds stage labels
        gene_col = "gene" if "gene" in df.columns else df.columns[1]

        union_genes = _select_top_by_stage(
            df, stage_col, gene_col, pval_cutoff, lfc_cutoff, top_n
        )
        gene_dict[sample] = set(union_genes)
        print(f"  {sample}: collected {len(union_genes)} genes")

    return gene_dict
# Step 2: 计算交集数据（逻辑不变）
def calculate_intersections(gene_dict, sample_names):
    intersections = {}
    selected_combinations = []  # 存储感兴趣的组合

    # 两两交集
    for comb in combinations(sample_names, 2):
        selected_combinations.append(comb)

    # 三个 MG 样本交集
    mg_samples = ["R-MG", "S-MG", "M-MG"]
    selected_combinations.append(tuple(mg_samples))

    # 所有六个样本的交集
    selected_combinations.append(tuple(sample_names))

    # 计算感兴趣的交集
    for comb in selected_combinations:
        intersect_set = set.intersection(*(gene_dict[sample] for sample in comb))
        intersections[comb] = len(intersect_set)

    return intersections

# Step 3: 绘制拼接图像（文件命名和标题做了调整）
def plot_combined(intersections, sample_names, save_path="trajectory_dependent_gene_intersection.png"):
    # 准备数据
    comb_keys = list(intersections.keys())
    counts = list(intersections.values())
    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    comb_keys = [comb_keys[i] for i in sorted_indices]
    counts = [counts[i] for i in sorted_indices]

    # 创建图像布局
    fig = plt.figure(figsize=(20, 10))  # 增加宽度以拉长间距
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[2, 1], hspace=0.05)

    # 上方柱状图
    ax_bar = fig.add_subplot(grid[0, 0])
    bar_positions = range(len(comb_keys))
    bar_width = 0.6

    # 使用 viridis 配色方案
    cmap = plt.cm.get_cmap("plasma", len(comb_keys))
    colors = [cmap(i) for i in range(len(comb_keys))]

    bars = ax_bar.bar(bar_positions, counts, width=bar_width, color=colors)
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylabel("Number of Genes")
    ax_bar.set_title("Intersection of Genes Along Pseudotime", fontsize=16)
    ax_bar.set_xlim(-0.5, len(comb_keys) - 0.5)  # 去掉左右空白

    # 在柱子上标注基因数量
    for bar, count in zip(bars, counts):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 5,
                    str(count), ha='center', fontsize=10)

    # 下方交集点线图
    ax_upset = fig.add_subplot(grid[1, 0])

    # 设置背景颜色交替
    for i in range(len(sample_names)):
        ax_upset.axhspan(i - 0.5, i + 0.5, color='lightgrey' if i % 2 == 0 else 'white', alpha=0.5)

    for i, comb in enumerate(comb_keys):
        for j, sample in enumerate(sample_names):
            if sample in comb:
                ax_upset.scatter(i, j, color='black', s=50)
        if len(comb) > 1:
            indices = [sample_names.index(sample) for sample in comb]
            ax_upset.plot([i, i], [min(indices), max(indices)], color='black', lw=1)

    ax_upset.set_xticks(bar_positions)  # 与柱状图共享X轴间距
    ax_upset.set_xticklabels([])
    ax_upset.set_yticks(range(len(sample_names)))
    ax_upset.set_yticklabels(sample_names)
    ax_upset.set_xlim(-0.5, len(comb_keys) - 0.5)  # 去掉左右空白
    ax_upset.set_ylim(-0.5, len(sample_names) - 0.5)

    # 保存图像（PNG和PDF）
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.savefig("trajectory_dependent_gene_intersection.pdf", dpi=400, bbox_inches='tight')
    plt.close()
    print(f"Combined plot saved as {save_path}")

# 主函数
def main(file_paths, sample_names):
    gene_dict = load_genelists(file_paths, sample_names)
    intersections = calculate_intersections(gene_dict, sample_names)
    plot_combined(intersections, sample_names, save_path="trajectory_dependent_gene_intersection.png")

# 运行主函数
if __name__ == "__main__":
    main(file_paths, sample_names)


# In[13]:


def calculate_and_export_intersections(gene_dict, sample_names, output_csv="trajectory_dependent_gene_intersections.csv"):
    intersections = {}

    # 计算两两交集
    for sample1, sample2 in combinations(sample_names, 2):
        intersect_set = gene_dict[sample1] & gene_dict[sample2]
        intersection_name = f"{sample1} & {sample2}"
        intersections[intersection_name] = list(intersect_set)

    # 计算 M-MG, R-MG, S-MG 三样本交集
    three_intersection_set = gene_dict["M-MG"] & gene_dict["R-MG"] & gene_dict["S-MG"]
    intersections["M-MG & R-MG & S-MG"] = list(three_intersection_set)

    # 计算所有六个样本的交集
    six_intersection_set = set.intersection(*(gene_dict[sample] for sample in sample_names))
    intersections["All Samples Intersection"] = list(six_intersection_set)

    # 转换为 DataFrame，并对齐不同长度的列表（不足部分以空字符串填充）
    max_length = max(len(genes) for genes in intersections.values())
    export_data = {name: genes + [""] * (max_length - len(genes)) for name, genes in intersections.items()}
    df = pd.DataFrame(export_data)
    df.to_csv(output_csv, index=False)
    print(f"Intersections saved to {output_csv}")

# 主函数
def main(file_paths, sample_names):
    gene_dict = load_genelists(file_paths, sample_names)
    calculate_and_export_intersections(gene_dict, sample_names, output_csv="trajectory_dependent_gene_intersections.csv")

if __name__ == "__main__":
    main(file_paths, sample_names)


# In[14]:


import pandas as pd

# 文件路径设置
novel_gene_file = "/data02/sunxuebo/project/scrnaseq/no-mammal/mammalnewgene/mammalian_new_genes_final.csv"  # 包含新基因的 CSV，列名为 Novel_Gene
intersections_file = "trajectory_dependent_gene_intersections.csv"
output_csv = "trajectory_dependent_novel_gene_intersections.csv"

# Step 1: 读取新基因 CSV，并提取 Novel_Gene 列构成集合
novel_df = pd.read_csv(novel_gene_file)
# 过滤空值并转换为字符串，构成新基因集合
novel_genes = set(novel_df['mouse'].dropna().astype(str).tolist())
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




