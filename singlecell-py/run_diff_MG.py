#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

# 预设全局变量，可以调整是否使用三物种共有HVGs
use_three_species_hvgs = True  # 设置为 True 使用三物种共有HVGs，False 使用样本间共有HVGs

##############################################################################
# 1. 读取三物种数据，并添加必要的 obs 信息
##############################################################################
path_m = "../../M-MG/1.subset/M-MG_cleaned.h5ad"
path_r = "../../R-MG/1.subset/R-MG_cleaned.h5ad"
path_s = "../../S-MG/1.subset/S-MG_cleaned.h5ad"

adata_m = sc.read_h5ad(path_m)
adata_r = sc.read_h5ad(path_r)
adata_s = sc.read_h5ad(path_s)

# 确保每个数据集都有 species 和 stage 信息
adata_m.obs['species'] = "M"
adata_r.obs['species'] = "R"
adata_s.obs['species'] = "S"
# 假设每个数据集的 obs['stage'] 已经包含 'stage1','stage2','stage3','stage4'

# 将三个数据集放入字典，方便后续调用
adata_dict = {"M": adata_m, "R": adata_r, "S": adata_s}

##############################################################################
# 2. 计算每个物种对在每个 stage 和 sample 的 Pearson 相关性
##############################################################################
species_pairs = [("M", "R"), ("M", "S"), ("R", "S")]
stages = ["stage1", "stage2", "stage3", "stage4"]

# 用于保存相关性数据
stage_corr_summary = []
sample_corr_summary = []

for sp1, sp2 in species_pairs:
    # 根据全局变量判断是否选择三物种共享的 HVGs
    if use_three_species_hvgs:
        hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
        hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
        hvgs_sp3 = adata_dict["S"].var_names[adata_dict["S"].var['highly_variable']]
        common_hvgs = list(set(hvgs_sp1).intersection(hvgs_sp2).intersection(hvgs_sp3))
        print(f"{sp1} vs {sp2} vs S 共有 HVGs 数量: {len(common_hvgs)}")
    else:
        hvgs_sp1 = adata_dict[sp1].var_names[adata_dict[sp1].var['highly_variable']]
        hvgs_sp2 = adata_dict[sp2].var_names[adata_dict[sp2].var['highly_variable']]
        common_hvgs = list(set(hvgs_sp1).intersection(set(hvgs_sp2)))
        print(f"{sp1} vs {sp2} 共有 HVGs 数量: {len(common_hvgs)}")

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
    
    # 计算每个 stage 的 Pearson 相关性
    stage_corr_vals = []
    for st in stages:
        mask_sp1 = (adata_pair.obs['species_batch'] == sp1) & (adata_pair.obs['stage'] == st)
        mask_sp2 = (adata_pair.obs['species_batch'] == sp2) & (adata_pair.obs['stage'] == st)
        adata_sp1 = adata_pair[mask_sp1]
        adata_sp2 = adata_pair[mask_sp2]
        
        if adata_sp1.n_obs == 0 or adata_sp2.n_obs == 0:
            stage_corr_vals.append(np.nan)
        else:
            mean_sp1 = np.asarray(adata_sp1.X.mean(axis=0)).flatten()
            mean_sp2 = np.asarray(adata_sp2.X.mean(axis=0)).flatten()
            mean_sp1 = np.nan_to_num(mean_sp1, nan=np.nanmean(mean_sp1))
            mean_sp2 = np.nan_to_num(mean_sp2, nan=np.nanmean(mean_sp2))
            corr_val, _ = pearsonr(mean_sp1, mean_sp2)
            stage_corr_vals.append(corr_val)
    
    # 构造阶段相关性的 DataFrame
    stage_corr_df = pd.DataFrame([stage_corr_vals], columns=stages)
    stage_corr_summary.append({
        'comparison': f"{sp1}_vs_{sp2}",
        'stage_correlation': stage_corr_df
    })

    # 计算每个 sample 的 Pearson 相关性
    sample_corr_vals = []
    for sample1 in adata_sp1.obs_names:
        for sample2 in adata_sp2.obs_names:
            sample_sp1 = adata_sp1[sample1]
            sample_sp2 = adata_sp2[sample2]
            mean_sp1 = np.asarray(sample_sp1.X.mean(axis=0)).flatten()
            mean_sp2 = np.asarray(sample_sp2.X.mean(axis=0)).flatten()
            corr_val, _ = pearsonr(mean_sp1, mean_sp2)
            sample_corr_vals.append(corr_val)
    
    # 构造样本相关性的 DataFrame
    sample_corr_df = pd.DataFrame([sample_corr_vals], columns=[f"{sp1}_samples", f"{sp2}_samples"])
    sample_corr_summary.append({
        'comparison': f"{sp1}_vs_{sp2}",
        'sample_correlation': sample_corr_df
    })

# 绘制 Stage Pearson 相关性图
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4), sharey=True)
for idx, (sp1, sp2) in enumerate(species_pairs):
    ax = axes[idx]
    corr_df = stage_corr_summary[idx]['stage_correlation']
    sns.heatmap(corr_df, cmap='RdBu_r', annot=True, fmt=".2f", vmin=-1, vmax=1, ax=ax)
    ax.set_title(f"{sp1} vs {sp2} - Stage Correlation", fontsize=12)
    ax.set_xticklabels(stages, rotation=45, ha='center')
    ax.set_yticks([])

fig.suptitle("Pearson Correlation Across Stages", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("Pearson_correlation_stages.png", dpi=300)
plt.show()

# 绘制 Sample Pearson 相关性图
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,4), sharey=True)
for idx, (sp1, sp2) in enumerate(species_pairs):
    ax = axes[idx]
    corr_df = sample_corr_summary[idx]['sample_correlation']
    sns.heatmap(corr_df, cmap='RdBu_r', annot=True, fmt=".2f", vmin=-1, vmax=1, ax=ax)
    ax.set_title(f"{sp1} vs {sp2} - Sample Correlation", fontsize=12)

fig.suptitle("Pearson Correlation Across Samples", fontsize=14)
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.savefig("Pearson_correlation_samples.png", dpi=300)
plt.show()
