#!/usr/bin/env python
# coding: utf-8
"""
Compute cosine similarity between cell‑type average expression profiles of two glands
(MG vs AG) **after gene‑wise (row‑wise) Z‑score standardisation on the *normalised* matrix**.

Changes (2025‑04‑24):
    • Removed `sc.pp.scale()` – single‑cell level scaling no longer applied.
    • Added *dotplot* visualisation of the cosine similarity matrix (user request).

Pipeline:
    1. Load AnnData → normalise counts → `log1p` → select HVGs.
    2. Build `newcelltype` (celltype‑gland) labels.
    3. Compute mean expression for each label → genes×groups matrix.
    4. Row‑wise Z‑score.
    5. Split columns into the two glands → pairwise cosine similarity (cross‑gland only).
    6. Save CSV, heatmap, **dotplot**.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cosine
from pathlib import Path

# -------------------- parameters --------------------
gland1 = "MG"     # first gland label in adata.obs['gland']
gland2 = "AG"     # second gland label

dataset      = "S-MAG"        # dataset name
celltype_col = "newcelltype"  # column with cell‑type (to be concatenated)
gland_col    = "gland"        # column with gland info
n_top_genes  = 5_000           # HVGs to keep
adata_path   = Path(f"../1.harmony/{dataset}_cleaned.h5ad")

# -------------------- load & basic pre‑proc --------------------
adata = sc.read_h5ad(adata_path)
adata.X = adata.layers["counts"]          # raw counts to X for normalisation
sc.pp.normalize_total(adata)              # CPM‑like scaling
sc.pp.log1p(adata)                        # log‑transform
sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
adata = adata[:, adata.var.highly_variable]
# (no sc.pp.scale call)

# create unique celltype‑gland labels
adata.obs[celltype_col] = (
    adata.obs[celltype_col].astype(str) + "-" + adata.obs[gland_col].astype(str)
)

# -------------------- build mean‑expression matrix --------------------
unique_groups = adata.obs[celltype_col].unique().tolist()
expr_mat = np.zeros((adata.n_vars, len(unique_groups)), dtype=np.float32)
for j, grp in enumerate(unique_groups):
    cells = adata.obs_names[adata.obs[celltype_col] == grp]
    expr_mat[:, j] = np.asarray(adata[cells].X.mean(axis=0)).ravel()
expr_df = pd.DataFrame(expr_mat, index=adata.var_names, columns=unique_groups)

# -------------------- row‑wise (gene) Z‑score --------------------
mu = expr_df.mean(axis=1)
sigma = expr_df.std(axis=1)
sigma[sigma == 0] = 1.0  # guard against division by zero
expr_z = expr_df.sub(mu, axis=0).div(sigma, axis=0)

# -------------------- split into glands --------------------
expr_z_g1 = expr_z[[c for c in expr_z.columns if c.endswith(gland1)]]
expr_z_g2 = expr_z[[c for c in expr_z.columns if c.endswith(gland2)]]

# -------------------- cosine similarity --------------------
results = []
for col_g1 in expr_z_g1.columns:
    for col_g2 in expr_z_g2.columns:
        vec1 = expr_z_g1[col_g1].values
        vec2 = expr_z_g2[col_g2].values
        cos_sim = 1 - cosine(vec1, vec2)
        results.append({
            f"cell_type_{gland1}": col_g1,
            f"cell_type_{gland2}": col_g2,
            "cosine_similarity": cos_sim,
        })

sim_df = pd.DataFrame(results)
outfile_csv = f"{dataset}_{gland1}_vs_{gland2}_cosine_similarity_zscore.csv"
sim_df.to_csv(outfile_csv, index=False)
print(f"[saved] {outfile_csv}")

# -------------------- heatmap --------------------
heat = sim_df.pivot(index=f"cell_type_{gland1}",
                    columns=f"cell_type_{gland2}",
                    values="cosine_similarity")
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(heat, aspect="auto")
ax.set_xticks(np.arange(len(heat.columns)))
ax.set_yticks(np.arange(len(heat.index)))
ax.set_xticklabels(heat.columns, rotation=90)
ax.set_yticklabels(heat.index)
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        val = heat.iat[i, j]
        if not np.isnan(val):
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6)
plt.colorbar(im)
plt.title(f"{dataset}: {gland1} vs {gland2} cell‑type cosine (row‑wise Z‑score)")
plt.tight_layout()
heat_out = f"{dataset}_{gland1}_vs_{gland2}_cosine_heatmap_zscore.png"
plt.savefig(heat_out, dpi=300)
print(f"[saved] {heat_out}")
plt.close()

# -------------------- dotplot --------------------
# replicate user‑provided styling, with safe scaling for point sizes
uniq_g1 = sorted(sim_df[f"cell_type_{gland1}"].unique())
uniq_g2 = sorted(sim_df[f"cell_type_{gland2}"].unique())

x_pos = sim_df[f"cell_type_{gland2}"].apply(lambda x: uniq_g2.index(x)).values
y_pos = sim_df[f"cell_type_{gland1}"].apply(lambda y: uniq_g1.index(y)).values
colors = sim_df["cosine_similarity"].values
# ensure sizes are positive: scale into 20‑120 range
sizes = 20 + 100 * ((colors - colors.min()) / (colors.max() - colors.min() + 1e-9))

fig, ax = plt.subplots(figsize=(8, 8))
scatter = ax.scatter(x_pos, y_pos, c=colors, s=sizes, cmap="jet", alpha=0.8)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label("Cosine Similarity", rotation=270, labelpad=10)

ax.set_xticks(range(len(uniq_g2)))
ax.set_yticks(range(len(uniq_g1)))
ax.set_xticklabels(uniq_g2, rotation=90)
ax.set_yticklabels(uniq_g1)
ax.set_xlabel(f"cell_type_{gland2}")
ax.set_ylabel(f"cell_type_{gland1}")
plt.title(f"Dotplot of Cosine Similarity in {dataset} ({gland1} vs {gland2})")
plt.tight_layout()
dot_out = f"{dataset}_{gland1}_vs_{gland2}_cosine_similarity_dotplot.png"
plt.savefig(dot_out, dpi=300)
print(f"[saved] {dot_out}")
plt.close()

