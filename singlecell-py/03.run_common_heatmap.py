#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Batch heat-map of Z-scored gene expression for a SINGLE AnnData file,
split by 'gland', using a celltype→gene dictionary (selection).

Author: <your name>
Date  : 2025-11-05
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# 1️⃣ 你的 celltype→Top5 基因字典
# (请将你上一步 'build_selection' 脚本输出的 'selection' 字典粘贴到这里)
# ------------------------------------------------------------------
selection = {
        "B lymphocytes": ["Blnk", "Tec", "Ralgps2", "Atp8a1", "Bmp2k"],
    "Basal epithelial cells": ["Krt5", "Krt17", "Col17a1", "Lamb3", "Trp63"],
    "Endothelial cells": ["Adgrl4", "Cdh5", "Flt1", "Myct1", "Mmrn2"],
    "Fibroblasts": ["Col1a1", "Dcn", "Col3a1", "Col5a1", "Mmp2"],
    "Innate immune cells": ["Csf1r", "Tyrobp", "C1qb", "C1qa", "Fcer1g"],
    "Keratinized epithelial cells": ["Rab25", "Sbsn", "Lypd3", "S100a14", "Krtdap"],
    "Luminal epithelial cells": ["Rasef", "Arhgef38", "Cldn3", "Cldn8", "Krt7"],
    "Proliferating epithelial cells": ["Ube2c", "Aurkb", "Top2a", "Ccnb1", "Cdk1"],
    "T lymphocytes": ["Itk", "Cd247", "Tox", "Cd28", "Icos"],
    "Vascular smooth muscle cells": ["Trpc6", "Gucy1b1", "Gucy1a1", "Aoc3", "Pde3a"]
}

# 衍生列表，保持顺序
celltype_list = list(selection.keys())
gene_list     = [g for ct in celltype_list for g in selection[ct]] 

print(f"[INFO] 成功加载 'selection' 字典。")
print(f"[INFO] 总共 {len(celltype_list)} 个细胞类型, {len(gene_list)} 个基因。")


# ------------------------------------------------------------------
# 2️⃣ 其他参数 (已修改)
# ------------------------------------------------------------------
adata_root = Path(r"./")
# ↓↓↓ 你的单个 .h5ad 文件名 (包含所有腺体)
H5AD_FILE  = "../../h5ad2seurat/processed/v9_all_log1p.h5ad" 

# ↓↓↓ 预处理后要绘图的腺体
GLANDS_TO_PLOT = ["EG", "MG", "SG"] 

# ↓↓↓ AnnData.obs 中用于细胞类型和腺体的列名
group_key  = "anno"          # (已修改) obs 中细胞类型列
gland_key  = "gland"         # obs 中腺体列

# (已恢复) 恢复你原来的配色方案
cmap       = sns.light_palette("#404040", n_colors=256, as_cmap=True) 
dpi        = 300

sns.set(style="white", font_scale=0.8)

# ------------------------------------------------------------------
# 3️⃣ 辅助函数 (用于预处理)
# ------------------------------------------------------------------
def rename_and_merge_glands(adata, gland_col):
    """
    应用: 1) SG -> EG; 2) AG -> SG, CG -> SG
    """
    if gland_col not in adata.obs.columns:
        raise ValueError(f"缺少 adata.obs['{gland_col}'] 列。")
    
    print(f"[INFO] 正在重命名/合并 '{gland_col}' 列...")
    # 确保是字符串，以便替换
    if pd.api.types.is_categorical_dtype(adata.obs[gland_col]):
        adata.obs[gland_col] = adata.obs[gland_col].astype(str)

    # Step 1: SG -> EG
    adata.obs[gland_col] = adata.obs[gland_col].replace({'SG': 'EG'})
    # Step 2: merge AG/CG into SG
    adata.obs[gland_col] = adata.obs[gland_col].replace({'AG': 'SG', 'CG': 'SG'})

    # 转回 Category
    adata.obs[gland_col] = adata.obs[gland_col].astype("category")
    print(f"[OK] 腺体合并完成。")
    print(f"[INFO] 重命名后 '{gland_col}' 计数：\n{adata.obs[gland_col].value_counts().to_string()}")
    return adata

# (MaSC 合并函数已移除)

# ------------------------------------------------------------------
# 4️⃣ 主流程 (已修改)
# ------------------------------------------------------------------

print("\n" + "="*30)
print("     开始执行热图绘制流程")
print("="*30)

# --- 1. 加载和预处理 (在循环外) ---
adata_path = adata_root / H5AD_FILE
if not adata_path.exists():
    raise FileNotFoundError(adata_path)

print(f"[INFO] 正在加载主 AnnData: {adata_path} ...")
adata_main = sc.read_h5ad(adata_path)
print(f"[OK] 主 AnnData 加载: {adata_main.n_obs:,} 个细胞, {adata_main.n_vars:,} 个基因。")

# 运行你的预处理
adata_main = rename_and_merge_glands(adata_main, gland_key)
# (MaSC 合并步骤已移除)

# 假设 .X 已经是归一化数据 (基于你的上一个脚本)
print(f"[INFO] 将使用 adata.X 作为归一化表达数据。")

# 提取颜色映射 (如果存在)
ct_color_map = None
if f"{group_key}_colors" in adata_main.uns:
    print("[INFO] 找到了细胞类型颜色映射。")
    # 确保使用 .obs 中的 categories 来正确匹配颜色
    if pd.api.types.is_categorical_dtype(adata_main.obs[group_key]):
        ct_color_map = dict(zip(
            adata_main.obs[group_key].cat.categories,
            adata_main.uns[f"{group_key}_colors"]
        ))
    else:
        print("[WARN] group_key 不是 category 类型，颜色映射可能不准确。")
        # 尝试基于 unique 值创建，但不保证顺序
        ct_color_map = dict(zip(
            np.unique(adata_main.obs[group_key]),
            adata_main.uns[f"{group_key}_colors"]
        ))


# --- 2. 主循环 (遍历腺体) ---
for gland in GLANDS_TO_PLOT:
    print(f"\n[INFO] Processing Gland: {gland} …")
    
    # (新增) 拆分 (subset)
    adata_gland = adata_main[adata_main.obs[gland_key] == gland].copy()
    
    if adata_gland.n_obs == 0:
        print(f"  [WARN] 在 {gland_key} 中找不到 {gland} 的细胞，跳过。")
        continue

    # ——— 基因 & 细胞类型实际可用列表 ————————————————
    genes_avail     = [g for g in gene_list if g in adata_gland.var_names]
    celltypes_avail = [ct for ct in celltype_list if ct in adata_gland.obs[group_key].unique()]
    
    if not genes_avail or not celltypes_avail:
        print(f"  [WARN] {gland} 中缺少足够的基因或细胞类型，跳过。")
        continue
    
    print(f"  [INFO] 找到 {len(celltypes_avail)} 个细胞类型和 {len(genes_avail)} 个基因。")

    # ——— 构建平均表达矩阵 gene × celltype ————————————
    expr = pd.DataFrame(index=genes_avail, columns=celltypes_avail, dtype=float)

    for ct in celltypes_avail:
        cells_mask = adata_gland.obs[group_key] == ct
        adata_subset = adata_gland[cells_mask, genes_avail]
        mean_vec = np.asarray(adata_subset.X.mean(axis=0)).ravel()
        expr.loc[genes_avail, ct] = mean_vec

    # ——— Z-score (行) & 旋转 ——————————
    expr_mean = expr.mean(axis=1)
    expr_std = expr.std(axis=1).replace(0, 1) # 替换0标准差
    
    expr_z = expr.sub(expr_mean, axis=0)
    expr_z = expr_z.div(expr_std, axis=0)
    
    heat   = expr_z.T  # 旋转: 行=celltype, 列=gene

    # 保证 row/col 顺序完全按 selection (的子集)
    heat = heat.loc[celltypes_avail, genes_avail]

    # ——— 行颜色条 ————————————————
    row_colors = None
    if ct_color_map:
        row_colors = [ct_color_map.get(ct, '#808080') for ct in heat.index] # 使用.get以防万一

    # ——— 画图 ——————————————
    figsize = (max(6, len(genes_avail) * 0.35),
               max(3, len(celltypes_avail) * 0.4))

    print(f"  [INFO] 正在绘制热图 (FigSize: {figsize[0]:.1f} x {figsize[1]:.1f})...")
    g = sns.clustermap(
        heat,
        row_cluster=False, col_cluster=False, 
        cmap=cmap,                          # (恢复) 灰色调
        vmin=0, vmax=2,                     # (恢复) 0-2 范围
        row_colors=row_colors,
        linewidths=0.1,
        figsize=figsize
    )

    g.ax_heatmap.set_xlabel("Gene")
    g.ax_heatmap.set_ylabel("Cell type")
    g.ax_heatmap.set_title(gland, pad=12) # 标题为腺体名
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

    out_name = f"{gland}_zscore_heatmap.pdf"
    g.savefig(out_name, dpi=dpi, bbox_inches="tight")
    plt.close()
    print(f"  [OK] Saved → {out_name}")

print("\n✅ All heatmaps generated!")
