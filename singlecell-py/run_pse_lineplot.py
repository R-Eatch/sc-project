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


# In[ ]:


# ========= 参数区域 =========
data_path   = "/data01/sunxuebo/project/scrnaseq/v8-python/"
datasetlist = ['M-MG', 'R-MG', 'S-MG', 'R-AG', 'R-CG', 'S-AG']
genelist    = ["Stc2",
                "Mboat1",
                "Gpx3",
                "Fkbp11",
                "Fam20a",
                "Sec11c",
                "Lipa",
                "Ssr4",
                "Plet1",
                "Slc28a3",
                "Arhgef38",
                "Cpeb4",
                "Arhgap26",
                "Rarres1",
                "Anxa1",
                "Pip",
                "Cidea",
               "Pigr",
               "Plin2",
               "Gata3",
               "Lgals7"
            ]         # ← 目标基因
n_stages    = 20                                     # pse_stage 个数
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

    # 仅对目标基因执行 scale，可进一步省内存
    sc.pp.scale(adata, zero_center=True, max_value=None)

    # ------ 划分 pse_stage ------
    adata.obs[group_key] = pd.qcut(
        adata.obs["pseudotime"],
        q=n_stages,
        labels=labels,
        duplicates="drop"
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
    expr = pd.DataFrame(
        adata[:, genelist].X,      # 这里的 X 已是 scaled
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
        estimator="median",
        errorbar=("pi", 90),        # 若不要阴影，改 None
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
plt.savefig("gene_pse_panel.png", dpi=300, bbox_inches="tight")
plt.show()


# In[ ]:




