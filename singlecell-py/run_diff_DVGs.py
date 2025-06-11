#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from itertools import combinations
import scanpy as sc
import matplotlib.pyplot as plt


# In[20]:


file_paths = [
    'D:/pythonfile/R-MG-LUMSEC.csv',
    'D:/pythonfile/R-AG-LUMSEC2.csv',
    'D:/pythonfile/S-MG-LUMSEC.csv',
    'D:/pythonfile/S-AG-LUMSEC2.csv',
    'D:/pythonfile/M-MG-LUMSEC1.csv',
    'D:/pythonfile/R-CG-LUMSEC2.csv'
]
sample_names = ["R-MG", "R-AG", "S-MG", "S-AG", "M-MG", "R-CG"]

def load_hvgs(file_paths, sample_names):
    hvgs_dict = {}
    for i, file_path in enumerate(file_paths):
        print(f"Loading sample: {sample_names[i]} from {file_path}...")
        df = pd.read_csv(file_paths[i])
        df1=df[df[df.columns[2]] < 0.01]
        df1[df.columns[0]].tolist()
        hvgs = df1[df.columns[0]].tolist() 
        hvgs_dict[sample_names[i]] = set(hvgs)  
        print(f"Sample {sample_names[i]}: {len(hvgs)} HVGs found.\n")
    return hvgs_dict
def calculate_intersections(hvgs_dict, sample_names):
    intersections = {}
    for r in range(1, len(sample_names) + 1):
        for combination in combinations(sample_names, r):
            intersect_set = set.intersection(*(hvgs_dict[sample] for sample in combination))
            intersections[combination] = len(intersect_set)
    return intersections

def plot_combined(intersections, sample_names, save_path="combined_plot_DVGs.png"):
    combinations = list(intersections.keys())
    counts = list(intersections.values())
    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    combinations = [combinations[i] for i in sorted_indices]
    counts = [counts[i] for i in sorted_indices]

    fig = plt.figure(figsize=(20, 10))  # 增加宽度以拉长间距
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[2, 1], hspace=0.05)

    ax_bar = fig.add_subplot(grid[0, 0])
    bar_positions = range(len(combinations))
    bar_width = 0.6
    ax_bar.bar(bar_positions, counts, width=bar_width, color='steelblue')  # 使用一致的 X 轴间距
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([])
    ax_bar.set_ylabel("Number of DVGs")
    ax_bar.set_title("Intersection of DVGs")
    ax_bar.set_xlim(-0.5, len(combinations) - 0.5)  # 去掉左右空白

    # 下方交集点线图
    ax_upset = fig.add_subplot(grid[1, 0])
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
    plt.close()
    print(f"Combined plot saved as {save_path}")
def calculate_and_export_intersections(hvgs_dict, sample_names, output_csv="DVGS-pairwise_intersections.csv"):
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
def main(file_paths, sample_names):
    hvgs_dict = load_hvgs(file_paths, sample_names)
    intersections = calculate_intersections(hvgs_dict, sample_names)
    plot_combined(intersections, sample_names)
    calculate_and_export_intersections(hvgs_dict, sample_names, output_csv="DVGs-pairwise_intersections.csv")

# 运行主函数
main(file_paths, sample_names)


# In[22]:


import os
import csv

###############################################################################
# 1) 定义每个样本对应的 gseGO, gseKEGG 结果文件
###############################################################################
# 注意：请根据自己的实际文件名进行修改
sample_files = {
    "M-MG": {
        "go":   "D:/111/DVGS-GSEA_Results/M-MG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/M-MG_gseKEGG_Result.csv"
    },
    "R-MG": {
        "go":   "D:/111/DVGS-GSEA_Results/R-MG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/R-MG_gseKEGG_Result.csv"
    },
    "S-MG": {
        "go":   "D:/111/DVGS-GSEA_Results/S-MG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/S-MG_gseKEGG_Result.csv"
    },
    "R-AG": {
        "go":   "D:/111/DVGS-GSEA_Results/R-AG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/R-AG_gseKEGG_Result.csv"
    },
    "S-AG": {
        "go":   "D:/111/DVGS-GSEA_Results/S-AG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/S-AG_gseKEGG_Result.csv"
    },
    "R-CG": {
        "go":   "D:/111/DVGS-GSEA_Results/R-CG_gseGO_Result.csv",
        "kegg": "D:/111/DVGS-GSEA_Results/R-CG_gseKEGG_Result.csv"
    },
}

# 所有样本名
all_samples = list(sample_files.keys())  # ["M-MG", "R-MG", "S-MG", "R-AG", "S-AG", "R-CG"]

###############################################################################
# 2) 读取 CSV 函数：返回一个 {ID: Description} 的字典
###############################################################################
def read_id_description(csv_path):
    """
    读取指定路径的 CSV 文件（包含列 "ID" 和 "Description"），
    返回一个字典 { id_value: description_value, ... }。
    若文件不存在或列缺失，返回空字典。
    """
    if not os.path.isfile(csv_path):
        print(f"[Warning] File not found: {csv_path}")
        return {}

    result_dict = {}
    with open(csv_path, "r", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        # 检查是否含有 ID, Description 两列
        if not {"ID", "Description"}.issubset(set(reader.fieldnames or [])):
            print(f"[Warning] File {csv_path} does not have 'ID'/'Description' columns.")
            return {}

        for row in reader:
            _id = row["ID"].strip()
            _desc = row["Description"].strip()
            if _id:
                result_dict[_id] = _desc
    return result_dict

###############################################################################
# 3) 一次性加载所有样本的 GO/KEGG 数据
###############################################################################
go_data_all   = {}  # go_data_all[sample_name] = {id: desc, ...}
kegg_data_all = {}  # kegg_data_all[sample_name] = {id: desc, ...}

for sample in all_samples:
    go_file   = sample_files[sample]["go"]
    kegg_file = sample_files[sample]["kegg"]
    go_data_all[sample]   = read_id_description(go_file)
    kegg_data_all[sample] = read_id_description(kegg_file)

###############################################################################
# 4) 取交集函数
###############################################################################
def intersect_id_description(dict_list):
    """
    给定若干个 {ID: Description} 字典，先求 ID 的交集，
    再返回一个同样形式的字典 {ID: Description}。
    Description 可以从第一个字典里拿或做更多处理，这里简化从第一个拿即可。
    """
    if not dict_list:
        return {}

    # 1) 先求所有字典的 ID 集合
    list_of_sets = [set(d.keys()) for d in dict_list if d]
    if not list_of_sets:
        return {}

    common_ids = set.intersection(*list_of_sets)  # 交集
    if not common_ids:
        return {}

    # 2) 从第一个字典中筛选出这些 ID，并保留其 Description
    #    如果需要更精细处理，可再检查每个字典里是否描述一致
    first_dict = dict_list[0]
    out_dict = {}
    for _id in common_ids:
        if _id in first_dict:
            out_dict[_id] = first_dict[_id]
    return out_dict

###############################################################################
# 5) 保存 {ID: Description} 到 CSV
###############################################################################
def save_id_desc_dict(data_dict, out_csv):
    """
    将 {ID: Description} 的字典写入 CSV，含表头 ID, Description。
    若字典为空，则打印提示，不写文件。
    """
    if not data_dict:
        print(f"[Info] Intersection/Difference is empty. Skip writing: {out_csv}")
        return

    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["ID", "Description"])
        for _id, _desc in data_dict.items():
            writer.writerow([_id, _desc])
    print(f"[OK] Wrote {len(data_dict)} rows to {out_csv}")

###############################################################################
# 6) 取差集函数：比较两个 {ID:Desc} 字典，返回只在第一个字典中的 ID
###############################################################################
def difference_id_description(dict_a, dict_b):
    """
    返回只在 dict_a 中、但不在 dict_b 中的 ID:Description。
    """
    out = {}
    if not dict_a:
        return out
    set_a = set(dict_a.keys())
    set_b = set(dict_b.keys())
    diff_ids = set_a - set_b
    for _id in diff_ids:
        out[_id] = dict_a[_id]
    return out

###############################################################################
# 7) 开始按照需求计算交集与差集
###############################################################################

# 7.1) 三样本交集：{M-MG, R-MG, S-MG}
tri_samples = ["M-MG", "R-MG", "S-MG"]
tri_go_dicts   = [go_data_all[s]   for s in tri_samples]
tri_kegg_dicts = [kegg_data_all[s] for s in tri_samples]

go_3 = intersect_id_description(tri_go_dicts)
kegg_3 = intersect_id_description(tri_kegg_dicts)

save_id_desc_dict(go_3,   "go_intersect_M-MG_R-MG_S-MG.csv")
save_id_desc_dict(kegg_3, "kegg_intersect_M-MG_R-MG_S-MG.csv")


# 7.2) 三组两两交集: (R-MG,R-AG), (R-MG,R-CG), (S-MG,S-AG)
pairs_for_intersect = [
    ("R-MG", "R-AG"),
    ("R-MG", "R-CG"),
    ("S-MG", "S-AG"),
]

for (s1, s2) in pairs_for_intersect:
    # GO intersection
    go_inter = intersect_id_description([go_data_all[s1], go_data_all[s2]])
    out_csv_go = f"go_intersect_{s1}_{s2}.csv"
    save_id_desc_dict(go_inter, out_csv_go)

    # KEGG intersection
    kegg_inter = intersect_id_description([kegg_data_all[s1], kegg_data_all[s2]])
    out_csv_kegg = f"kegg_intersect_{s1}_{s2}.csv"
    save_id_desc_dict(kegg_inter, out_csv_kegg)

# 7.3) 六样本交集: {M-MG, R-MG, S-MG, R-AG, S-AG, R-CG}
all_go_dicts   = [go_data_all[s]   for s in all_samples]
all_kegg_dicts = [kegg_data_all[s] for s in all_samples]

go_6 = intersect_id_description(all_go_dicts)
kegg_6 = intersect_id_description(all_kegg_dicts)

save_id_desc_dict(go_6,   "go_intersect_all6.csv")
save_id_desc_dict(kegg_6, "kegg_intersect_all6.csv")

###############################################################################
# 8) 额外需求：计算三组 (S-MG,S-AG), (R-MG,R-AG), (R-MG,R-CG) 
#    在 KEGG 中各自独有的通路（差集）
###############################################################################
pairs_for_diff = [
    ("S-MG", "S-AG"),
    ("R-MG", "R-AG"),
    ("R-MG", "R-CG"),
]

for (s1, s2) in pairs_for_diff:
    dict1 = kegg_data_all[s1]
    dict2 = kegg_data_all[s2]

    # s1 独有 (in s1 but not in s2)
    diff_s1 = difference_id_description(dict1, dict2)
    csv_s1 = f"kegg_unique_{s1}_vs_{s2}.csv"
    save_id_desc_dict(diff_s1, csv_s1)

    # s2 独有 (in s2 but not in s1)
    diff_s2 = difference_id_description(dict2, dict1)
    csv_s2 = f"kegg_unique_{s2}_vs_{s1}.csv"
    save_id_desc_dict(diff_s2, csv_s2)

###############################################################################
print("[Done] All intersections and differences have been processed.")


# In[ ]:




