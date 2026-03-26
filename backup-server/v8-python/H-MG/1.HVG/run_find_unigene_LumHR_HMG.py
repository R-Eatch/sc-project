#!/usr/bin/env python
# coding: utf-8
"""
find_LumHR_unique_genes_v3.py
───────────────────────────────────
• 计算并保存三物种 LumHR DEG（完整 & 过滤后）
• 统计三物种所有细胞的基因表达百分比
• 获取人乳腺 LumHR DEG 的 human‑unique 基因：
    - logFC > LOGFC_MIN, p_adj < P_ADJ_MAX     （LumHR‑specific in human）
    - 该基因在 M/R/S 任意一个物种的所有细胞中均 ≤5% 细胞表达
      （即所有物种 >95% 细胞不表达）
• 对三个 LumHR 亚群分别输出结果

⚠️ 与 v2 的主要区别
──────────────────
去除了「三物种 LumHR‑specific 基因合集」功能（即不再汇总 `union_deg`）。
"""

from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# ------------ 输入路径 ------------
PATHS = {
    "M": "/data01/sunxuebo/project/scrnaseq/v8-python/M-MG/1.subset/M-MG_cleaned.h5ad",
    "R": "/data01/sunxuebo/project/scrnaseq/v8-python/R-MG/1.subset/R-MG_cleaned.h5ad",
    "S": "/data01/sunxuebo/project/scrnaseq/v8-python/S-MG/1.subset/S-MG_cleaned.h5ad",
}
HUMAN_DEG_CSV = "H-MG_ranked_genes_author_cell_type-mouse.csv"
HR_GROUPS = ["LummHR-SCGB", "LummHR-active", "LummHR-major"]

# ------------ 阈值参数 ------------
LOGFC_MIN    = 0.25   # 最小 logFC
P_ADJ_MAX    = 0.05   # 最大调整后 P 值
PCT_HR_MIN   = 0.10   # LumHR 内至少 10% 细胞表达
PCT_REST_MAX = 0.01   # 非 LumHR ≤5% 细胞表达
EXP_PCT_MAX  = 0.05   # 全细胞表达百分比阈值（5%）

# ------------ 输出目录 ------------
OUT = Path("LumHR_unique_results")
OUT.mkdir(exist_ok=True)

# --------------------------------------------------
# ❶ 三物种：计算 LumHR DEG & 统计全细胞表达百分比
# --------------------------------------------------
union_var   = set()                  # 三物种全部基因合集
gene_pcts   = defaultdict(list)      # {gene: [pct_M, pct_R, pct_S]}

for sp, fp in PATHS.items():
    print(f"[{sp}] 读取 {fp}")
    ad = sc.read_h5ad(fp)
    union_var.update(ad.var_names)

    # ---------- ① LumHR vs rest DEG ----------
    sc.tl.rank_genes_groups(
        ad,
        groupby="newcelltype",
        groups=["LumHR"],
        reference="rest",
        layer="normalized" if "normalized" in ad.layers else None,
        use_raw=False,
        method="wilcoxon",
        pts=True,
    )

    de = sc.get.rank_genes_groups_df(ad, group="LumHR").rename(columns={
        "names":            "gene",
        "logfoldchanges":   "logFC",
        "pct_nz_group":     "pct_HR",
        "pct_nz_reference": "pct_rest",
    })

    # 保存完整表
    de.to_csv(OUT / f"{sp}_LumHR_DEG_full.csv", index=False)

    keep = de.query(
        "logFC > @LOGFC_MIN and pvals_adj < @P_ADJ_MAX "
        "and pct_HR >= @PCT_HR_MIN and pct_rest <= @PCT_REST_MAX"
    )
    keep.to_csv(OUT / f"{sp}_LumHR_DEG_filtered.csv", index=False)
    print(f"  ↳ LumHR‑specific genes (filtered): {len(keep)}")

    # ---------- ② 全细胞表达百分比 ----------
    X = ad.layers.get("counts", ad.X)  # 计数矩阵
    if sparse.issparse(X):
        n_cells_expr = (X > 0).sum(axis=0).A1
    else:
        n_cells_expr = (X > 0).sum(axis=0)

    pct_expr = n_cells_expr / ad.n_obs
    for gene, pct in zip(ad.var_names, pct_expr):
        gene_pcts[gene].append(pct)

print(f"\n✓ Union var_names across M/R/S: {len(union_var)}")

# ---------- ③ 计算三物种无表达基因 ----------
nonexpr_genes = {
    g for g, pcts in gene_pcts.items()
    if len(pcts) == len(PATHS) and all(p < EXP_PCT_MAX for p in pcts)
}
print(f"✓ Genes unexpressed (≤5% cells) in ALL species: {len(nonexpr_genes)}\n")

# --------------------------------------------------
# ❷ 人乳腺 LumHR DEG → human‑unique & non‑expressed in M/R/S
# --------------------------------------------------
human = pd.read_csv(HUMAN_DEG_CSV)

for hr in HR_GROUPS:
    df_hr = human[human["group"] == hr].copy()

    # 阈值过滤
    df_hr = df_hr.query(
        "logfoldchanges > @LOGFC_MIN and pvals_adj < @P_ADJ_MAX"
    )

    # 保留三物种均无表达的基因
    uniq_raw = df_hr[df_hr["mouse_gene"].isin(nonexpr_genes)]
    uniq_raw.to_csv(OUT / f"{hr}_human_unique_nonexpr_raw.csv", index=False)

    # 与 union_var 交集（确保 downstream 可用）
    uniq_inData = uniq_raw[uniq_raw["mouse_gene"].isin(union_var)]
    uniq_inData.to_csv(OUT / f"{hr}_human_unique_nonexpr_inData.csv", index=False)

    print(
        f"{hr}: after DEG filter {len(df_hr)} genes"
        f" → unique_nonexpr_raw {len(uniq_raw)}, inData {len(uniq_inData)}"
    )

print("\n所有结果已写入", OUT.resolve())
