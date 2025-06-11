#!/usr/bin/env python
# coding: utf-8

# In[2]:


# mg_apocrine_gene_discovery.py (patched 2025‑05‑30)
"""
Identify mammary‑gland‑specific genes (MG) absent in all apocrine glands (AG + CG → "AP")
across three species (M, R, S).

⚠️ Patch notes
──────────────
* **Fix Scanpy ValueError** when a species lacks AP cells (e.g. mouse).  
  ‑ If a species has no AP category, we now fall back to a simple **detection‑rate filter**
  ‑ Treat its AP set as empty (i.e. no apocrine expression measured).
* Logic & output files unchanged.
"""

# ------------------------------------------------------------
# ❶ Imports & global parameters
# ------------------------------------------------------------
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# — Hard‑coded cleaned h5ad paths —
PATHS = {
    "M-MG": "/data01/sunxuebo/project/scrnaseq/v8-python/M-MG/1.subset/M-MG_cleaned.h5ad",
    "R-MG": "/data01/sunxuebo/project/scrnaseq/v8-python/R-MG/1.subset/R-MG_cleaned.h5ad",
    "S-MG": "/data01/sunxuebo/project/scrnaseq/v8-python/S-MG/1.subset/S-MG_cleaned.h5ad",
    "R-AG": "/data01/sunxuebo/project/scrnaseq/v8-python/R-AG/1.subset/R-AG_cleaned.h5ad",
    "S-AG": "/data01/sunxuebo/project/scrnaseq/v8-python/S-AG/1.subset/S-AG_cleaned.h5ad",
    "R-CG": "/data01/sunxuebo/project/scrnaseq/v8-python/R-CG/1.subset/R-CG_cleaned.h5ad",
}

SPECIES = ["M", "R", "S"]
GLANDS  = ["MG", "AG", "CG"]

# Thresholds (tweak freely)
LOGFC_MIN  = 0.25
P_ADJ_MAX  = 0.05
PCT_MG_MIN = 0.05
PCT_AP_MAX = 0.05

# — Output dir —
OUT_DIR = Path("./results")
OUT_DIR.mkdir(exist_ok=True)


# In[ ]:


# ------------------------------------------------------------
# ❷ Read data & concatenate
# ------------------------------------------------------------
if not os.path.exists('all_epi.h5ad'):
    print("Loading h5ad files …")
    adata_list, keys = [], []
    for tag, p in PATHS.items():
        path = Path(p)
        if not path.exists():
            print("  ⚠️  Missing", path)
            continue
        ad = sc.read_h5ad(path)
        sp, gl = tag.split("-")
        ad.obs["species"] = sp
        ad.obs["gland"]   = gl
        ad.layers.setdefault("normalized", ad.X.copy())
        adata_list.append(ad)
        keys.append(tag)
        print("  ✓", tag, ad.n_obs, "cells")
    
    adata_all = sc.concat(adata_list, join="outer", label="samplebatch", keys=keys)
    adata_all.obs["gland3"] = adata_all.obs["gland"].replace({"AG":"AP", "CG":"AP", "MG":"MG"})
    adata_all.write('all_epi.h5ad')
else:
    print("已存在 all_epi.h5ad，skip generate")
    adata_all=sc.read_h5ad('all_epi.h5ad')
    print("Loading h5ad files …")


# In[ ]:


# ------------------------------------------------------------
# ❸ MG vs AP gene sets per species
# ------------------------------------------------------------
print("Running per‑species MG/AP parsing …")
print("Running per-species MG/AP parsing …")
mg_sets, ap_sets = {}, {}
for sp in SPECIES:
    ad_sp = adata_all[adata_all.obs["species"] == sp].copy()
    ad_sp.raw = None  # 防止 layer+raw 冲突
    cats = ad_sp.obs["gland3"].unique().tolist()

    # ---------- 情况 A: 该物种同时存在 MG 和 AP ----------
    if {"MG", "AP"}.issubset(cats):
        sc.tl.rank_genes_groups(
            ad_sp,
            groupby="gland3",
            groups=["MG"],
            reference="AP",
            layer="normalized",
            use_raw=False,
            method="wilcoxon",
            pts=True,
        )
        de = sc.get.rank_genes_groups_df(ad_sp, group="MG")

        # --- 1) 统一列名 ---
        if "pct_nz_group" in de.columns:  # Scanpy ≥1.9
            de = de.rename(columns={
                "pct_nz_group": "pct_MG",
                "pct_nz_reference": "pct_AP",
                "logfoldchanges": "logfc",
            })
            de.to_csv(f'{sp}_gland_DEG.csv')
        elif "pct_nz" in de.columns:      # 旧版本 (≤1.8) 没有 reference 列
            # 在这种版本里无法直接得 pct_AP，需手动算
            de = de.rename(columns={
                "pct_nz": "pct_MG",
                "logfoldchanges": "logfc",
            })

        # --- 2) 若仍缺检测率列，手动计算 ---
        if "pct_AP" not in de.columns:
            pct_ref = (ad_sp[ad_sp.obs["gland3"] == "AP"].layers["normalized"] > 0).mean(axis=0).A1
            de["pct_AP"] = pd.Series(pct_ref, index=ad_sp.var_names).reindex(de["names"]).values
        if "pct_MG" not in de.columns:
            pct_grp = (ad_sp[ad_sp.obs["gland3"] == "MG"].layers["normalized"] > 0).mean(axis=0).A1
            de["pct_MG"] = pd.Series(pct_grp, index=ad_sp.var_names).reindex(de["names"]).values

        # --- 3) 过滤得到集合 ---
        mg_sets[sp] = set(
            de.query(
                "logfc > @LOGFC_MIN & pvals_adj < @P_ADJ_MAX & "
                "pct_MG >= @PCT_MG_MIN & pct_AP <= @PCT_AP_MAX"
            )["names"]
        )
        ap_sets[sp] = set(de.query("pct_AP > @PCT_AP_MAX")["names"])
        print(
            f"  {sp}: Wilcoxon – {len(mg_sets[sp])} MG-specific, {len(ap_sets[sp])} AP-expressed"
        )

    # ---------- 情况 B: 该物种没有 AP 细胞 (如 mouse) ----------
    else:
        mg_cells = ad_sp[ad_sp.obs["gland3"] == "MG"]
        nz = (mg_cells.layers["normalized"] > 0)  # 稀疏布尔矩阵
        pct_mg = nz.mean(axis=0).A1
        pct_mg = pd.Series(pct_mg, index=mg_cells.var_names)
        mg_sets[sp] = set(pct_mg[pct_mg >= PCT_MG_MIN].index)
        ap_sets[sp] = set()  # 无 AP 信息
        print(f"  {sp}: no AP cells – detection filter → {len(mg_sets[sp])} MG-expressed genes")


# In[ ]:


# ------------------------------------------------------------
# ④ Presence / absence matrix + core & pairwise gene sets
# ------------------------------------------------------------
mg_cols = [f"{sp}_MG" for sp in SPECIES]
ap_cols = [f"{sp}_AP" for sp in SPECIES]

all_genes = sorted(adata_all.var_names)
mat = pd.DataFrame(False, index=all_genes, columns=mg_cols + ap_cols)
for sp in SPECIES:
    mat.loc[list(mg_sets[sp]), f"{sp}_MG"] = True      # 乳腺
    mat.loc[list(ap_sets[sp]), f"{sp}_AP"] = True      # 顶泌腺

mat.to_csv(OUT_DIR / "mg_presence_absence_matrix.csv")

core_mask = mat[mg_cols].all(axis=1) & (~mat[ap_cols].any(axis=1))
core_genes = mat.index[core_mask].tolist()

pairwise_genes = set()
for sp1, sp2 in [("M","R"), ("M","S"), ("R","S")]:
    mask = mat[[f"{sp1}_MG", f"{sp2}_MG"]].all(axis=1) & (~mat[ap_cols].any(axis=1))
    pairwise_genes.update(mat.index[mask])
pairwise_genes = sorted(pairwise_genes - set(core_genes))

pd.Series(core_genes, name="gene").to_csv(OUT_DIR/"mg_core_MRS.csv", index=False)
pd.Series(pairwise_genes, name="gene").to_csv(OUT_DIR/"mg_pairwise.csv", index=False)
print(f"→ CORE (M-R-S): {len(core_genes)} genes → mg_core_MRS.csv")
print(f"→ PAIRWISE (any two): {len(pairwise_genes)} genes → mg_pairwise.csv")

