#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# -*- coding: utf-8 -*-
"""
单图展示多基因在 pse_stage 中的表达变化
"""

import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings, re, gc
from pathlib import Path
import math
import numpy as np


# In[ ]:


# ========= 参数区域 =========
data_path   = "/data01/sunxuebo/project/scrnaseq/v8-python/"
datasetlist = ['M-MG', 'R-MG', 'S-MG', 'R-AG', 'R-CG', 'S-AG']
genelist    = ["Stc2",
                "Atp1b1",
		"Pdk4",
		"Gata3"
            ]         # ← 目标基因
n_stages    = 15                                     # pse_stage 个数
group_key   = "pse_stage"
celltype_key = "newcelltype"                        # 细胞类型列名

# 针对每个数据集的细胞类型筛选；空列表 = 不筛选
celltype_filter = {
    'M-MG': [],
    'R-MG': [],
    'S-MG': [],
    'R-AG': [],
    'R-CG': [],
    'S-AG': []
}
# 若希望用统一颜色（如示例中的物种色），可在此自定义
dataset_palette = sns.color_palette("Set2", n_colors=len(datasetlist))
# ============================


# In[ ]:





# In[ ]:


sns.set(style="ticks", font_scale=1.1)
labels   = [f"pse_{i+1}" for i in range(n_stages)]
df_all   = []     # 汇总表


for ds, color in zip(datasetlist, dataset_palette):
    ad_path = Path(data_path) / ds / "9.stage_DEG" / f"{ds}_cleaned_pse.h5ad"
    print(f"── 处理：{ds}")
    adata = sc.read_h5ad(ad_path)

    # ------ 细胞类型筛选 ------
    ct_list = celltype_filter.get(ds, [])
    if ct_list:     # 非空才做筛选
        mask = adata.obs[celltype_key].isin(ct_list)
        adata = adata[mask].copy()
        print(f"   · 细胞类型筛选后：{adata.n_obs} cells")

    if adata.n_obs == 0:
        warnings.warn(f"{ds} 无可用细胞，跳过")
        continue

    # ------ 归一化层 → adata.X ------
    if "normalized" in adata.layers:
        adata.X = adata.layers["normalized"].copy()
    else:
        warnings.warn(f"{ds} 未找到 normalized 层，直接用 adata.X")

    # NOTE: 不再对全体基因 sc.pp.scale(adata)。
    # 我们将在抽取 genelist 表达后，仅对这些目标基因做 gene-wise z-score（Scaled expression）。

    # ------ 划分 pse_stage ------
    # 1) pseudotime 在每个数据集内归一化到 0-1
    if "pseudotime" not in adata.obs.columns:
        warnings.warn(f"{ds} 未找到 obs['pseudotime']，跳过")
        continue

    pt = pd.to_numeric(adata.obs["pseudotime"], errors="coerce").astype(float)
    good = np.isfinite(pt.values)
    if good.sum() == 0:
        warnings.warn(f"{ds} pseudotime 全为 NA/inf，跳过")
        continue
    if good.sum() != adata.n_obs:
        adata = adata[good].copy()
        pt = pt[good]

    pt_min, pt_max = float(pt.min()), float(pt.max())
    if pt_max == pt_min:
        warnings.warn(f"{ds} pseudotime 无跨度（min==max），跳过")
        continue

    pt01 = (pt - pt_min) / (pt_max - pt_min)
    adata.obs["pseudotime_01"] = pt01.values

    # 2) 在 [0,1] 上固定等宽切成 20 个 bins（pse_1 ... pse_20）
    edges = np.linspace(0.0, 1.0, n_stages + 1)
    adata.obs[group_key] = pd.cut(
        adata.obs["pseudotime_01"],
        bins=edges,
        labels=labels,
        include_lowest=True
    )

    stages_sorted = sorted(
        adata.obs[group_key].unique(),
        key=lambda x: int(re.sub(r"[^0-9]", "", str(x)))
    )
    adata.obs[group_key] = pd.Categorical(
        adata.obs[group_key],
        categories=stages_sorted,
        ordered=True
    )

    # ------ 抽取表达矩阵到 DataFrame ------
    # 只对目标基因做 scaled expression（每基因在该数据集内 across cells 的 z-score）
    Xg = adata[:, genelist].X
    # 稀疏矩阵转 dense，便于逐基因计算均值/方差
    if hasattr(Xg, "toarray"):
        Xg = Xg.toarray()
    else:
        Xg = np.asarray(Xg)

    mu = Xg.mean(axis=0, keepdims=True)
    sd = Xg.std(axis=0, keepdims=True)
    sd[sd == 0] = 1.0
    Xg_z = (Xg - mu) / sd

    expr = pd.DataFrame(
        Xg_z,
        columns=genelist,
        index=adata.obs_names
    )
    expr[group_key] = adata.obs[group_key].values
    expr["sample"]  = ds
    df_all.append(expr)

    # 立即释放内存
    del adata
    gc.collect()

# ------------------ 汇总绘图 ------------------


# In[ ]:





# In[ ]:


if not df_all:
    raise RuntimeError("所有数据集中都未找到目标基因，无法绘图")

df_all = pd.concat(df_all, axis=0)


# In[ ]:


# ---------- 筛选要画的基因 ----------
plot_genes = genelist               # 你已决定直接用原列表
n_genes = len(plot_genes)

# ---------- 子图排版 ----------
n_cols = 4
n_rows = math.ceil(n_genes / n_cols)
fig, axes = plt.subplots(n_rows, n_cols,
                         figsize=(3*n_cols, 3*n_rows),
                         sharey=False)

axes = axes.flatten()               # 保证可迭代

for ax, gene in zip(axes, plot_genes):
    sns.lineplot(
        data=df_all,
        x=group_key,
        y=gene,
        hue="sample",
        estimator="mean",
        errorbar=None,              # 不画误差阴影
        palette=dataset_palette,
        linewidth=1,                # 细线
        marker="o",
        markersize=3,               # 小点
        ax=ax
    )
    ax.set_title(gene, style='italic', fontsize=10)
    ax.set_xlabel("Pseudotime bins")
    ax.set_ylabel("Scaled expression")
    ax.tick_params(axis='x', rotation=0, length=2)  # 小刻度
    ax.set_xticklabels([])          # 删掉 pse_1 … 文本
    sns.despine(ax=ax)

# 若子图不满最后一行，关掉多余框
for extra_ax in axes[n_genes:]:
    extra_ax.remove()

# ---------- 合并图例 ----------
handles, labels_ = axes[0].get_legend_handles_labels()
fig.legend(handles, labels_, title="Sample",
           bbox_to_anchor=(1.02, 0.5), loc="center left")
for ax in axes[:n_genes]:
    ax.get_legend().remove()

fig.tight_layout()
plt.savefig("gene_pse_panel.pdf", dpi=300, bbox_inches="tight")
plt.show()


# In[ ]:





