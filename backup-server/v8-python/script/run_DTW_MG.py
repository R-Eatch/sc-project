#!/usr/bin/env python
# coding: utf-8

# In[1]:


#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dtw import dtw  # 使用dtw库进行动态时间规整

##############################################################################
# 新增四个开关变量，用于控制是否进行MG / AG 的 sample / stage 分析
##############################################################################
do_MG_sample = True   # 是否进行MG数据的sample分析
do_MG_stage  = True   # 是否进行MG数据的stage分析
do_AG_sample = True   # 是否进行AG数据的sample分析
do_AG_stage  = True   # 是否进行AG数据的stage分析

# 预设全局变量，可以调整是否使用三物种共有HVGs
use_three_species_hvgs = True  # 设置为 True 使用三物种共有HVGs，False 使用样本间共有HVGs

##############################################################################
# 1. 读取三物种数据，并添加必要的 obs 信息 (MG部分)
##############################################################################
path_m = "../../M-MG/1.subset/M-MG_cleaned.h5ad"
path_r = "../../R-MG/1.subset/R-MG_cleaned.h5ad"
path_s = "../../S-MG/1.subset/S-MG_cleaned.h5ad"

adata_m = sc.read_h5ad(path_m)
adata_r = sc.read_h5ad(path_r)
adata_s = sc.read_h5ad(path_s)

# 为每个数据集添加物种标签
adata_m.obs['species'] = "M"
adata_r.obs['species'] = "R"
adata_s.obs['species'] = "S"

# 将三个数据集放入字典，方便后续调用
adata_dict = {"M": adata_m, "R": adata_r, "S": adata_s}

##############################################################################
# 预设三个物种样本的时间顺序（手动定义）
##############################################################################
R_order = ['R-MG-E17', 'R-MG-E23','R-MG-P1' ,'R-MG-8WK-1',
  'R-MG-23WK-3', 'R-MG-23WK-4','R-MG-GES12', 'R-MG-GES17', 'R-MG-GES23','R-MG-LA','R-MG-LA-2']
S_order = ['S-MG-P7',
           'S-MG-P20', 
           'S-MG-3MTH',
           'S-MG-8M-3',
           'S-MG-8M-4',
           'S-MG-GES14',
           'S-MG-LA-1',
            'S-MG-LA-2']
M_order = ['M-MG-E13_5',
           'M-MG-E16_5', 'M-MG-P1','M-MG-3WK-1',
           'M-MG-3WK-2',
           'M-MG-8WK-1',
           'M-MG-8WK-2',
           'M-MG-GES13_5',
           'M-MG-GES16_5',
           'M-MG-LA-1',
           'M-MG-LA-2']

##############################################################################
# 2 & 3. 针对各物种对，进行“样本”级别的Manhattan距离 & DTW分析
#        以及绘制曼哈顿距离热图 (样本维度)
##############################################################################
if do_MG_sample:
    species_pairs = [("M", "R"), ("M", "S"), ("R", "S")]

    sample_distance_summary = {}
    dtw_summary = {}

    for sp1, sp2 in species_pairs:
        print(f"正在处理物种对: {sp1} vs {sp2}", flush=True)
        
        # 根据全局变量判断是否选择三物种共享的 HVGs
        if use_three_species_hvgs:
            hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
            hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
            hvgs_sp3 = adata_dict["S"].var_names[adata_dict["S"].var['highly_variable']]
            common_hvgs = list(set(hvgs_sp1).intersection(hvgs_sp2).intersection(hvgs_sp3))
            print(f"{sp1} vs {sp2} vs S 共有 HVGs 数量: {len(common_hvgs)}", flush=True)
        else:
            hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
            hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
            common_hvgs = list(set(hvgs_sp1).intersection(set(hvgs_sp2)))
            print(f"{sp1} vs {sp2} 共有 HVGs 数量: {len(common_hvgs)}", flush=True)
        
        # 针对当前 pair，从原始数据中只保留共同 HVGs
        a1 = adata_dict[sp1][:, common_hvgs].copy()
        a2 = adata_dict[sp2][:, common_hvgs].copy()
        
        # 合并两个数据集，利用 concatenate 保证批次区分
        adata_pair = a1.concatenate(
            a2,
            batch_key="species_batch",
            batch_categories=[sp1, sp2],
            join='inner'
        )
        
        # 根据预设的样本顺序，确定各物种样本的顺序
        if sp1 == "R":
            sp1_samples = R_order
        elif sp1 == "S":
            sp1_samples = S_order
        elif sp1 == "M":
            sp1_samples = M_order
        else:
            sp1_samples = adata_pair.obs.loc[adata_pair.obs['species_batch'] == sp1, 'sample'].unique().tolist()
            
        if sp2 == "R":
            sp2_samples = R_order
        elif sp2 == "S":
            sp2_samples = S_order
        elif sp2 == "M":
            sp2_samples = M_order
        else:
            sp2_samples = adata_pair.obs.loc[adata_pair.obs['species_batch'] == sp2, 'sample'].unique().tolist()
        
        print(f"使用的 {sp1} 样本顺序: {sp1_samples}", flush=True)
        print(f"使用的 {sp2} 样本顺序: {sp2_samples}", flush=True)
        
        # 计算每个样本的均值表达向量
        sp1_sample_means = {}
        for sample in sp1_samples:
            mask = (adata_pair.obs['species_batch'] == sp1) & (adata_pair.obs['sample'] == sample)
            sp1_sample_means[sample] = np.asarray(adata_pair[mask].X.mean(axis=0)).flatten()
        
        sp2_sample_means = {}
        for sample in sp2_samples:
            mask = (adata_pair.obs['species_batch'] == sp2) & (adata_pair.obs['sample'] == sample)
            sp2_sample_means[sample] = np.asarray(adata_pair[mask].X.mean(axis=0)).flatten()
        
        # 计算样本间的Manhattan距离（绝对差值之和）
        sample_distance_mat = pd.DataFrame(index=sp1_samples, columns=sp2_samples, dtype=float)
        for sample1 in sp1_samples:
            for sample2 in sp2_samples:
                print(f"正在计算样本: {sample1} vs {sample2}", flush=True)
                vec1 = sp1_sample_means[sample1]
                vec2 = sp2_sample_means[sample2]
                manhattan_dist = np.sum(np.abs(vec1 - vec2))
                sample_distance_mat.loc[sample1, sample2] = manhattan_dist
        
        # 保存Manhattan距离矩阵到 CSV
        out_csv = f"Manhattan_distance_{sp1}_vs_{sp2}_samples.csv"
        sample_distance_mat.to_csv(out_csv)
        print(f"保存 Manhattan 距离矩阵到 {out_csv}", flush=True)
        sample_distance_summary[f"{sp1}_vs_{sp2}"] = sample_distance_mat.copy()
        
        # -------------------------------
        # DTW 分析：基于样本间 Manhattan 距离矩阵进行 DTW 对齐
        # -------------------------------
        seq_x = np.arange(len(sp1_samples))
        seq_y = np.arange(len(sp2_samples))
        
        def cost_function(x, y):
            i = int(x)
            j = int(y)
            return sample_distance_mat.iloc[i, j]
        
        dtw_distance, dtw_cost_matrix, dtw_acc_cost_matrix, dtw_path = dtw(seq_x, seq_y, dist=cost_function)
        path_x, path_y = dtw_path
        print(f"DTW 距离 ({sp1} vs {sp2}):", dtw_distance)
        
        dtw_summary[f"{sp1}_vs_{sp2}"] = {
            "dtw_distance": dtw_distance,
            "path": (path_x, path_y),
            "acc_cost_matrix": dtw_acc_cost_matrix
        }
        
        # 绘制DTW累积代价矩阵及最佳对齐路径 (在这里将"累计代价"改为"Accumulated Cost")
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(dtw_acc_cost_matrix.T, cmap='viridis',
                         cbar=True, cbar_kws={'label': 'Accumulated Cost'},
                         xticklabels=sp1_samples, yticklabels=sp2_samples)
        plt.xlabel(f"{sp1} Samples (ordered)")
        plt.ylabel(f"{sp2} Samples (ordered)")
        plt.title(f"DTW Alignment ({sp1} vs {sp2}): Accumulated Cost Matrix")
        
        # 在heatmap中绘制箭头，表示最优对齐路径 (optimal alignment path)
        for k in range(len(path_x) - 1):
            x0, y0 = path_x[k], path_y[k]
            dx = path_x[k+1] - x0
            dy = path_y[k+1] - y0
            plt.arrow(x0 + 0.5, y0 + 0.5, dx, dy, color='black', 
                      head_width=0.3, head_length=0.3, length_includes_head=True)
        
        plt.tight_layout()
        out_fig = f"DTW_alignment_heatmap_{sp1}_vs_{sp2}.png"
        plt.savefig(out_fig, dpi=300)
        plt.show()
        print(f"保存 DTW 对齐热图到 {out_fig}", flush=True)

    # 绘制“样本”级曼哈顿距离热图 (在这里将"曼哈顿距离"改为"Manhattan Distance")
    order_dict = {"R": R_order, "S": S_order, "M": M_order}
    for sp1, sp2 in species_pairs:
        sp1_samples = order_dict.get(sp1, adata_dict[sp1].obs['sample'].unique().tolist())
        sp2_samples = order_dict.get(sp2, adata_dict[sp2].obs['sample'].unique().tolist())
        sample_distance_mat = sample_distance_summary[f"{sp1}_vs_{sp2}"]
        
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(sample_distance_mat, cmap='plasma',
                         cbar_kws={'label': 'Manhattan Distance'},
                         xticklabels=sp2_samples, yticklabels=sp1_samples)
        plt.xlabel(f"{sp2} Samples")
        plt.ylabel(f"{sp1} Samples")
        plt.title(f"Manhattan Distance Heatmap (samples): {sp1} vs {sp2}")
        plt.tight_layout()
        out_fig = f"Manhattan_distance_heatmap_{sp1}_vs_{sp2}_samples.png"
        plt.savefig(out_fig, dpi=300)
        plt.show()
        print(f"保存曼哈顿距离热图 (样本)到 {out_fig}", flush=True)


##############################################################################
# 4. 基于 stage 分组进行分析（MG部分）
##############################################################################
if do_MG_stage:
    species_pairs = [("M", "R"), ("M", "S"), ("R", "S")]

    for sp1, sp2 in species_pairs:
        print(f"正在处理物种对 (stage分组): {sp1} vs {sp2}", flush=True)
        
        if use_three_species_hvgs:
            hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
            hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
            hvgs_sp3 = adata_dict["S"].var_names[adata_dict["S"].var['highly_variable']]
            common_hvgs = list(set(hvgs_sp1).intersection(hvgs_sp2).intersection(hvgs_sp3))
            print(f"{sp1} vs {sp2} vs S 共有 HVGs 数量: {len(common_hvgs)}", flush=True)
        else:
            hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
            hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
            common_hvgs = list(set(hvgs_sp1).intersection(set(hvgs_sp2)))
            print(f"{sp1} vs {sp2} 共有 HVGs 数量: {len(common_hvgs)}", flush=True)
        
        a1 = adata_dict[sp1][:, common_hvgs].copy()
        a2 = adata_dict[sp2][:, common_hvgs].copy()
        
        adata_pair = a1.concatenate(
            a2,
            batch_key="species_batch",
            batch_categories=[sp1, sp2],
            join='inner'
        )

        desired_stage_order = ["stage0", "stage1", "stage2", "stage3", "stage4"]
        
        sp1_all_stages = adata_pair.obs.loc[adata_pair.obs['species_batch'] == sp1, 'stage'].unique().tolist()
        sp1_stages = [stage for stage in desired_stage_order if stage in sp1_all_stages]
        
        sp2_all_stages = adata_pair.obs.loc[adata_pair.obs['species_batch'] == sp2, 'stage'].unique().tolist()
        sp2_stages = [stage for stage in desired_stage_order if stage in sp2_all_stages]

        sp1_stage_means = {}
        for stage in sp1_stages:
            mask = (adata_pair.obs['species_batch'] == sp1) & (adata_pair.obs['stage'] == stage)
            sp1_stage_means[stage] = np.asarray(adata_pair[mask].X.mean(axis=0)).flatten()
        
        sp2_stage_means = {}
        for stage in sp2_stages:
            mask = (adata_pair.obs['species_batch'] == sp2) & (adata_pair.obs['stage'] == stage)
            sp2_stage_means[stage] = np.asarray(adata_pair[mask].X.mean(axis=0)).flatten()
        
        stage_distance_mat = pd.DataFrame(index=sp1_stages, columns=sp2_stages, dtype=float)
        for stage1 in sp1_stages:
            for stage2 in sp2_stages:
                vec1 = sp1_stage_means[stage1]
                vec2 = sp2_stage_means[stage2]
                manhattan_dist = np.sum(np.abs(vec1 - vec2))
                stage_distance_mat.loc[stage1, stage2] = manhattan_dist
        
        out_csv_stage = f"Manhattan_distance_{sp1}_vs_{sp2}_stages.csv"
        stage_distance_mat.to_csv(out_csv_stage)
        print(f"保存曼哈顿距离矩阵 (stage) 到 {out_csv_stage}", flush=True)
        
        # DTW 分析 (stage)
        seq_x = np.arange(len(sp1_stages))
        seq_y = np.arange(len(sp2_stages))
        
        def cost_function_stage(x, y):
            i = int(x)
            j = int(y)
            return stage_distance_mat.iloc[i, j]
        
        dtw_distance, dtw_cost_matrix, dtw_acc_cost_matrix, dtw_path = dtw(seq_x, seq_y, dist=cost_function_stage)
        path_x, path_y = dtw_path
        print(f"DTW 距离 (stage) ({sp1} vs {sp2}):", dtw_distance)
        
        # 绘制 Accumulated Cost 矩阵 (stage)
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(dtw_acc_cost_matrix.T, cmap='viridis',
                         cbar_kws={'label': 'Accumulated Cost'},
                         xticklabels=sp1_stages, yticklabels=sp2_stages)
        plt.xlabel(f"{sp1} Stages (ordered)")
        plt.ylabel(f"{sp2} Stages (ordered)")
        plt.title(f"DTW Alignment (stage) ({sp1} vs {sp2}): Accumulated Cost Matrix")
        for k in range(len(path_x) - 1):
            x0, y0 = path_x[k], path_y[k]
            dx = path_x[k+1] - x0
            dy = path_y[k+1] - y0
            plt.arrow(x0 + 0.5, y0 + 0.5, dx, dy, color='black', 
                      head_width=0.3, head_length=0.3, length_includes_head=True)
        plt.tight_layout()
        out_fig_stage = f"DTW_alignment_heatmap_{sp1}_vs_{sp2}_stages.png"
        plt.savefig(out_fig_stage, dpi=300)
        plt.show()
        print(f"保存 DTW 对齐热图 (stage) 到 {out_fig_stage}", flush=True)
        
        # 额外绘制曼哈顿距离热图 (stage) -> "Manhattan Distance"
        plt.figure(figsize=(10, 8))
        ax = sns.heatmap(stage_distance_mat, cmap='plasma',
                         cbar_kws={'label': 'Manhattan Distance'},
                         xticklabels=sp2_stages, yticklabels=sp1_stages)
        plt.xlabel(f"{sp2} Stages")
        plt.ylabel(f"{sp1} Stages")
        plt.title(f"Manhattan Distance Heatmap (stage): {sp1} vs {sp2}")
        plt.tight_layout()
        out_fig_stage_manhattan = f"Manhattan_distance_heatmap_{sp1}_vs_{sp2}_stages.png"
        plt.savefig(out_fig_stage_manhattan, dpi=300)
        plt.show()
        print(f"保存曼哈顿距离热图 (stage) 到 {out_fig_stage_manhattan}", flush=True)


##############################################################################
# 5. 单独重新加载 R-AG 和 S-AG 数据，并对其进行 sample 和 stage 分组分析（AG部分）
##############################################################################
path_r_ag = "../../R-AG/1.subset/R-AG_cleaned.h5ad"
path_s_ag = "../../S-AG/1.subset/S-AG_cleaned.h5ad"

adata_rag = sc.read_h5ad(path_r_ag)
adata_sag = sc.read_h5ad(path_s_ag)

# 添加物种标签
adata_rag.obs['species'] = "R-AG"
adata_sag.obs['species'] = "S-AG"

adata_dict_ag = {"R-AG": adata_rag, "S-AG": adata_sag}

for grouping in ["sample", "stage"]:
    if grouping == "sample" and not do_AG_sample:
        continue
    if grouping == "stage" and not do_AG_stage:
        continue

    print(f"正在处理 R-AG vs S-AG, 分组依据: {grouping}", flush=True)
    
    hvgs_r = adata_rag.var_names[adata_rag.var['highly_variable']]
    hvgs_s = adata_sag.var_names[adata_sag.var['highly_variable']]
    common_hvgs = list(set(hvgs_r).intersection(set(hvgs_s)))
    print(f"R-AG vs S-AG 共有 HVGs 数量: {len(common_hvgs)}", flush=True)
    
    a_r = adata_rag[:, common_hvgs].copy()
    a_s = adata_sag[:, common_hvgs].copy()
    
    adata_pair_ag = a_r.concatenate(
         a_s,
         batch_key="species_batch",
         batch_categories=["R-AG", "S-AG"],
         join='inner'
    )
    
    if grouping == "sample":
        #groups_r = sorted(adata_pair_ag.obs.loc[adata_pair_ag.obs['species_batch'] == "R-AG", 'sample'].unique().tolist())
        groups_r = ['R-AG-E26', 'R-AG-P1''R-AG-8WK-1','R-AG-10WK-1','R-AG-25WK-1', 'R-AG-25WK-2','R-AG-GES12','R-AG-GES17','R-AG-GES23','R-AG-LA','R-AG-LA-2'
]
        #groups_s = sorted(adata_pair_ag.obs.loc[adata_pair_ag.obs['species_batch'] == "S-AG", 'sample'].unique().tolist())
        groups_s=['S-AG-P20','S-AG-3MTH-1','S-AG-3MTH-2','S-AG-8M-2','S-AG-8M-3','S-AG-GES14','S-AG-LA-1','S-AG-LA-2']
        group_label = "Samples"
    else:
        groups_r = sorted(adata_pair_ag.obs.loc[adata_pair_ag.obs['species_batch'] == "R-AG", 'stage'].unique().tolist())
        groups_s = sorted(adata_pair_ag.obs.loc[adata_pair_ag.obs['species_batch'] == "S-AG", 'stage'].unique().tolist())
        group_label = "Stages"
    
    print(f"使用的 R-AG {group_label} 顺序: {groups_r}", flush=True)
    print(f"使用的 S-AG {group_label} 顺序: {groups_s}", flush=True)
    
    group_means_r = {}
    for grp in groups_r:
        mask = (adata_pair_ag.obs['species_batch'] == "R-AG") & (adata_pair_ag.obs[grouping] == grp)
        group_means_r[grp] = np.asarray(adata_pair_ag[mask].X.mean(axis=0)).flatten()
    
    group_means_s = {}
    for grp in groups_s:
        mask = (adata_pair_ag.obs['species_batch'] == "S-AG") & (adata_pair_ag.obs[grouping] == grp)
        group_means_s[grp] = np.asarray(adata_pair_ag[mask].X.mean(axis=0)).flatten()
    
    distance_mat_ag = pd.DataFrame(index=groups_r, columns=groups_s, dtype=float)
    for grp_r in groups_r:
        for grp_s in groups_s:
            vec_r = group_means_r[grp_r]
            vec_s = group_means_s[grp_s]
            manhattan_dist = np.sum(np.abs(vec_r - vec_s))
            distance_mat_ag.loc[grp_r, grp_s] = manhattan_dist
    
    out_csv_ag = f"Manhattan_distance_R-AG_vs_S-AG_{grouping}.csv"
    distance_mat_ag.to_csv(out_csv_ag)
    print(f"保存曼哈顿距离矩阵 (R-AG vs S-AG, {grouping}) 到 {out_csv_ag}", flush=True)
    
    # DTW 分析
    seq_x = np.arange(len(groups_r))
    seq_y = np.arange(len(groups_s))
    
    def cost_func_ag(x, y):
        i = int(x)
        j = int(y)
        return distance_mat_ag.iloc[i, j]
    
    dtw_distance_ag, dtw_cost_mat_ag, dtw_acc_cost_mat_ag, dtw_path_ag = dtw(seq_x, seq_y, dist=cost_func_ag)
    path_x_ag, path_y_ag = dtw_path_ag
    print(f"DTW 距离 (R-AG vs S-AG, {grouping}):", dtw_distance_ag)
    
    # 绘制 Accumulated Cost 矩阵 (AG)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(dtw_acc_cost_mat_ag.T, cmap='viridis',
                     cbar_kws={'label': 'Accumulated Cost'},
                     xticklabels=groups_r, yticklabels=groups_s)
    plt.xlabel(f"R-AG {group_label} (ordered)")
    plt.ylabel(f"S-AG {group_label} (ordered)")
    plt.title(f"DTW Alignment (R-AG vs S-AG, {grouping}): Accumulated Cost Matrix")
    for k in range(len(path_x_ag) - 1):
        x0, y0 = path_x_ag[k], path_y_ag[k]
        dx = path_x_ag[k+1] - x0
        dy = path_y_ag[k+1] - y0
        plt.arrow(x0+0.5, y0+0.5, dx, dy, color='black', head_width=0.3, head_length=0.3, length_includes_head=True)
    plt.tight_layout()
    out_fig_ag = f"DTW_alignment_heatmap_R-AG_vs_S-AG_{grouping}.png"
    plt.savefig(out_fig_ag, dpi=300)
    plt.show()
    print(f"保存 DTW 对齐热图 (R-AG vs S-AG, {grouping}) 到 {out_fig_ag}", flush=True)
    
    # 额外绘制Manhattan Distance热图 (AG)
    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(distance_mat_ag, cmap='plasma',
                     cbar_kws={'label': 'Manhattan Distance'},
                     xticklabels=groups_s, yticklabels=groups_r)
    plt.xlabel(f"S-AG {group_label}")
    plt.ylabel(f"R-AG {group_label}")
    plt.title(f"Manhattan Distance Heatmap (R-AG vs S-AG, {grouping})")
    plt.tight_layout()
    out_fig_ag_manhattan = f"Manhattan_distance_heatmap_R-AG_vs_S-AG_{grouping}.png"
    plt.savefig(out_fig_ag_manhattan, dpi=300)
    plt.show()
    print(f"保存曼哈顿距离热图 (R-AG vs S-AG, {grouping}) 到 {out_fig_ag_manhattan}", flush=True)

print("所有计算完成！", flush=True)


# In[ ]:





# In[ ]:




