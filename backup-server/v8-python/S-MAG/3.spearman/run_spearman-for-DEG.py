#!/usr/bin/env python
# coding: utf-8

# In[ ]:

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import scipy
import seaborn as sns

gland1 = 'MG'
gland2 = 'AG'
dataset = 'S-MAG'
celltype_column = 'newcelltype'
gland_column = 'gland'
n_top_genes=2000
path1=f"/data01/sunxuebo/project/scrnaseq/v8-python/S-MG/1.subset/S-MG_cleaned-holo.h5ad"
path2=f"/data01/sunxuebo/project/scrnaseq/v8-python/S-AG/1.subset/S-AG_cleaned.h5ad"
vmin, vmax = 0.5,0.9
# In[4]:


adata1 = sc.read_h5ad(path1)
adata2 = sc.read_h5ad(path2)
adata1.obs['clusters']=adata1.obs['leiden']
adata2.obs['clusters']=adata2.obs['leiden']
adata1.obs = adata1.obs.drop(columns=["cellid"])
adata2.obs = adata2.obs.drop(columns=["cellid"])
adata = adata1.concatenate(adata2,batch_key = 'glandbatch')

adata.X = adata.layers['counts'].copy()

adata.obs['newcelltype'] = (adata.obs['newcelltype'].astype(str) + '-' +  adata.obs[gland_column].astype(str))

# In[ ]:


sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
adata = adata[:, adata.var.highly_variable]
#sc.pp.scale(adata, max_value=10)


# In[6]:


adata_A = adata[adata.obs[gland_column]==gland1]
adata_B = adata[adata.obs[gland_column]==gland2]



cell_types_A = adata_A.obs[celltype_column].unique()
cell_types_B = adata_B.obs[celltype_column].unique()


# In[9]:


correlation_results = []

for celltype_A in cell_types_A:
    for celltype_B in cell_types_B:
        if celltype_A == celltype_B:
            continue

        print(f"current celltype: {celltype_A} vs {celltype_B}")

        # 选出 gland1 + celltype_A
        subset_gland1 = adata[
            (adata.obs[celltype_column] == celltype_A) &
            (adata.obs[gland_column] == gland1)
        ]

        # 选出 gland2 + celltype_B
        subset_gland2 = adata[
            (adata.obs[celltype_column] == celltype_B) &
            (adata.obs[gland_column] == gland2)
        ]

        # 确保子集非空
        if subset_gland1.shape[0] > 0 and subset_gland2.shape[0] > 0:

            subset_gland1.X = np.nan_to_num(subset_gland1.X, nan=0.0)
            subset_gland2.X = np.nan_to_num(subset_gland2.X, nan=0.0)

            mean_expression_gland1 = np.asarray(subset_gland1.X.mean(axis=0)).flatten()
            mean_expression_gland2 = np.asarray(subset_gland2.X.mean(axis=0)).flatten()

            corr, pval = spearmanr(mean_expression_gland1, mean_expression_gland2,nan_policy='omit')
            correlation_results.append({
                'cell_type_A': celltype_A,
                'cell_type_B': celltype_B,
                'correlation': corr,
                'p-value': pval
            })

# 整理结果并过滤 p < 0.05
correlation_df = pd.DataFrame(correlation_results).dropna()
#correlation_df = correlation_df[correlation_df['p-value'] < 0.1]
print(correlation_df)


# In[ ]


# In[11]:


# In[2]:


#correlation_df=pd.read_csv("./S-MAG-pse_correlation_results.csv")


# In[10]:


# ===== 绘制热图部分 =====
# 1) 将长表格 (long format) 转成矩阵形式 (wide format)
#    行为 cell_type_A，列为 cell_type_B，矩阵值为 correlation
correlation_matrix = correlation_df.pivot(
    index='cell_type_A',
    columns='cell_type_B',
    values='correlation'
)

import re

def _numeric_stage_sort(labels):
    """
    Sort labels like 'pse_1-MG', 'pse_12-AG' … by the numeric part between '_' and '-'.
    """
    def _key(s):
        m = re.search(r'_(\d+)(?:-|$)', s)
        return int(m.group(1)) if m else float('inf')
    return sorted(labels, key=_key)
    
# 1) 构造行列全集并按数字排序
stage_x = _numeric_stage_sort(correlation_df['cell_type_B'].unique())
stage_y = _numeric_stage_sort(correlation_df['cell_type_A'].unique())

# 2) 生成全矩阵并保持 NaN（方便 imshow 上写值时跳过）
corr_matrix = (correlation_df
               .pivot(index='cell_type_A', columns='cell_type_B', values='correlation')
               .reindex(index=stage_y, columns=stage_x))

# 3) 动态尺寸 & 字体
fig_w_hm = max(4, len(stage_x) * 0.35 + 2)
fig_h_hm = max(4, len(stage_y) * 0.35 + 2)
tick_font = max(6, 12 - 0.25 * max(len(stage_x), len(stage_y)))
title_font = max(10, 16 - 0.30 * max(len(stage_x), len(stage_y)))

# 4) 开始绘制
fig, ax = plt.subplots(figsize=(fig_w_hm, fig_h_hm))
im = ax.imshow(corr_matrix, cmap='RdBu_r')
colors = correlation_df['correlation'].values
im.set_clim(colors.min(), colors.max())
ax.invert_yaxis()
# 坐标轴
ax.set_xticks(np.arange(len(stage_x)))
ax.set_yticks(np.arange(len(stage_y)))
ax.set_xticklabels(stage_x, rotation=90, ha='right', rotation_mode='anchor', fontsize=tick_font)
ax.set_yticklabels(stage_y, fontsize=tick_font)

# 在每格写数值（跳过 NaN）
for i in range(len(stage_y)):
    for j in range(len(stage_x)):
        val = corr_matrix.iat[i, j]
        if not np.isnan(val):
            ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                    color='black', fontsize=tick_font-1)

# 颜色条
cbar = fig.colorbar(im)
cbar.ax.tick_params(labelsize=tick_font)
cbar.set_label('Spearman ρ', rotation=270, labelpad=12, fontsize=tick_font)

# 标题与布局
ax.set_xlabel('cell_type_B', fontsize=tick_font)
ax.set_ylabel('cell_type_A', fontsize=tick_font)
plt.title(f'Spearman Correlation Heatmap in {dataset}', fontsize=title_font)
plt.tight_layout()
plt.savefig(f'{dataset}-DEG_correlation_heatmap.png', dpi=300)
plt.show()
# ===== 3) 若需要同时保存结果DataFrame =====
correlation_df.to_csv(f"{dataset}-DEG_correlation_results.csv", index=False)


# In[4]:


#correlation_df=pd.read_csv("./S-MAG-celltype_correlation_results.csv")
#correlation_df


# In[17]:


# 假设在脚本其他位置，你已获得 correlation_df，并且 import 了 matplotlib 等库


# In[5]:


# 1) 提取并按数字排序 pse‑stage（cell_type_A / cell_type_B）
unique_celltype_A = _numeric_stage_sort(correlation_df['cell_type_A'].unique())
unique_celltype_B = _numeric_stage_sort(correlation_df['cell_type_B'].unique())

# 2) 映射到数值索引（横纵坐标）
x_positions = correlation_df['cell_type_B'].apply(lambda x: unique_celltype_B.index(x))
y_positions = correlation_df['cell_type_A'].apply(lambda y: unique_celltype_A.index(y))

# ------------------------------------------------------------------
# 根据 stage 数量动态计算画布 & 字体大小
# ------------------------------------------------------------------
max_stage = max(len(unique_celltype_A), len(unique_celltype_B))

#   尺寸：每个 stage 给 ~0.4 inch，外加 2 inch 边框，至少 4×4 inch
fig_w = max(4, len(unique_celltype_B) * 0.4 + 2)
fig_h = max(4, len(unique_celltype_A) * 0.4 + 2)

#   字体：stage 越多越缩小；最小 6 pt
tick_font  = max(6, 12 - 0.25 * max_stage)
title_font = max(10, 16 - 0.30 * max_stage)

# 3) 颜色 & 点大小
colors = correlation_df['correlation'].values
sizes  = 20 * -np.log10(correlation_df['p-value'] + 1e-16)


# In[ ]:





# In[6]:


# ------------------------------------------------------------------
# Dot‑plot
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(fig_w, fig_h))

# Scatter plot – give the limits directly …
sc = ax.scatter(
    x_positions, y_positions,
    c=colors, s=sizes,
    cmap='jet', alpha=0.8
)
sc.set_clim(vmin, vmax)

# 右侧留出空间放 colorbar & p‑value legend
plt.subplots_adjust(right=0.78)

# colorbar
cbar = plt.colorbar(sc, ax=ax, fraction=0.08, pad=0.02)
cbar.set_label('Spearman ρ', rotation=270, labelpad=12, fontsize=tick_font)

# 坐标轴
ax.set_xticks(range(len(unique_celltype_B)))
ax.set_yticks(range(len(unique_celltype_A)))
ax.set_xticklabels(unique_celltype_B, rotation=90, fontsize=tick_font)
ax.set_yticklabels(unique_celltype_A, fontsize=tick_font)
ax.set_xlabel('cell_type_B', fontsize=tick_font)
ax.set_ylabel('cell_type_A', fontsize=tick_font)

# p‑value 大小图例
pvals = [0.05, 0.01, 0.001]
handles = [ax.scatter([], [], s=20 * -np.log10(p+1e-16),
                      color='gray', alpha=0.8) for p in pvals]
labels  = [f'{p:g}' for p in pvals]

ax.legend(
    handles, labels, title='p‑value',
    loc='lower right', bbox_to_anchor=(1.46, 0),
    frameon=False, fontsize=tick_font, title_fontsize=tick_font
)

ax.set_title(f'{dataset} | dot‑plot of stage correlations',
             fontsize=title_font, pad=20)
plt.tight_layout()
plt.savefig(f'{dataset}-DEG_correlation_dotplot.png', dpi=300)
plt.savefig(f'{dataset}-DEG_correlation_dotplot.pdf', dpi=300)
plt.close()


# In[ ]:





# In[ ]:





# In[ ]:




