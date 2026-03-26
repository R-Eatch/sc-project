#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Pairwise stage‑level correlation analysis for three scRNA‑seq datasets
────────────────────────────────────────────────────────────────────────
* Datasets:  M‑MG  |  R‑MG  |  S‑MG
* Metric:    Spearman correlation on mean expression of HVGs (n_top_genes = 2 000)
* Stage col: adata.obs['stage']
* Output:    • <ds1>_vs_<ds2>_stage_correlation.csv
             • <ds1>_vs_<ds2>_heatmap.png
             • <ds1>_vs_<ds2>_dotplot.png
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from itertools import combinations
import pathlib
import sys

# ────── User‑configurable section ────────────────────────────────────────────
DATASETS      = ["M-MG", "R-MG", "S-MG"]          # exactly three dataset IDs
CELL_COL      = "stage"                             # stage column in .obs
N_TOP_GENES   = 1000                                # HVG count
DATA_TEMPLATE = "../../{ds}/1.subset/{ds}_cleaned.h5ad"  # path template
OUT_DIR       = pathlib.Path("pairwise_stage_correlation")
# ─────────────────────────────────────────────────────────────────────────────

OUT_DIR.mkdir(parents=True, exist_ok=True)
print("Datasets →", DATASETS)
print("Output   →", OUT_DIR.absolute())

# ────── Helper functions ────────────────────────────────────────────────────

def load_dataset(ds_id: str) -> sc.AnnData:
    """Load a cleaned AnnData object with raw counts in .layers['counts']."""
    path = pathlib.Path(DATA_TEMPLATE.format(ds=ds_id))
    if not path.exists():
        sys.exit(f"✗ cannot find {path}")
    ad = sc.read_h5ad(path)
    if "cellid" in ad.obs.columns:
        ad.obs.drop(columns=["cellid"], inplace=True)
    return ad


def compute_pair(ds1: str, ds2: str):
    print(f"\n▶ analysing {ds1} ↔ {ds2}")
    ad1, ad2 = load_dataset(ds1), load_dataset(ds2)

    # concatenate to ensure identical gene space and joint HVG selection
    ad = ad1.concatenate(ad2, batch_key="dataset", batch_categories=[ds1, ds2])
    ad.X = ad.layers["counts"].copy()

    # preprocessing
    ad.obs[CELL_COL] = ad.obs[CELL_COL].astype(str)
    sc.pp.normalize_total(ad)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, n_top_genes=N_TOP_GENES)
    ad = ad[:, ad.var.highly_variable]

    # split back
    ad_A = ad[ad.obs["dataset"] == ds1]
    ad_B = ad[ad.obs["dataset"] == ds2]

    stages_A = ad_A.obs[CELL_COL].unique()
    stages_B = ad_B.obs[CELL_COL].unique()

    rows = []
    for sA in stages_A:
        for sB in stages_B:
            sub_A = ad[(ad.obs[CELL_COL] == sA) & (ad.obs["dataset"] == ds1)]
            sub_B = ad[(ad.obs[CELL_COL] == sB) & (ad.obs["dataset"] == ds2)]
            if sub_A.n_obs == 0 or sub_B.n_obs == 0:
                continue
            mA = np.asarray(sub_A.X.mean(axis=0)).flatten()
            mB = np.asarray(sub_B.X.mean(axis=0)).flatten()
            corr, pval = spearmanr(mA, mB, nan_policy="omit")
            rows.append({"cell_type_A": sA, "cell_type_B": sB, "correlation": corr, "p-value": pval})

    corr_df = pd.DataFrame(rows).dropna()
    tag = f"{ds1}_vs_{ds2}"
    corr_df.to_csv(OUT_DIR / f"{tag}_stage_correlation.csv", index=False)

    # heatmap
    mat = corr_df.pivot(index="cell_type_A", columns="cell_type_B", values="correlation")
    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(mat, cmap="coolwarm", vmin=0.4, vmax=0.7)
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=90)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index)
    for i in range(len(mat.index)):
        for j in range(len(mat.columns)):
            val = mat.iloc[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6)
    plt.colorbar(im, ax=ax)
    plt.title(f"Spearman Correlation – {tag}")
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"{tag}_heatmap.png", dpi=300)
    plt.close()

    # dot‑plot
    uniq_A = sorted(corr_df["cell_type_A"].unique())
    uniq_B = sorted(corr_df["cell_type_B"].unique())
    xpos = corr_df["cell_type_B"].apply(lambda x: uniq_B.index(x))
    ypos = corr_df["cell_type_A"].apply(lambda y: uniq_A.index(y))
    colors = corr_df["correlation"].values
    sizes  = 20 * -np.log10(corr_df["p-value"] + 1e-16)

    fig, ax = plt.subplots(figsize=(8, 8))
    scat = ax.scatter(xpos, ypos, c=colors, s=sizes, cmap="jet", alpha=0.8)
    scat.set_clim(colors.min(), colors.max())
    plt.colorbar(scat, ax=ax, fraction=0.3, pad=0.02, shrink=0.5).set_label("Spearman r", rotation=270, labelpad=10)
    for p in [0.05, 0.01, 0.001]:
        ax.scatter([], [], s=20 * -np.log10(p + 1e-16), color="gray", alpha=0.8, label=f"p ≤ {p}")
    ax.legend(title="p-value", loc="lower right", bbox_to_anchor=(1.35, 0), frameon=False)
    ax.set_xticks(range(len(uniq_B)))
    ax.set_xticklabels(uniq_B, rotation=90)
    ax.set_yticks(range(len(uniq_A)))
    ax.set_yticklabels(uniq_A)
    ax.set_xlabel(f"stage – {ds2}")
    ax.set_ylabel(f"stage – {ds1}")
    plt.title(f"Dot‑plot of stage correlation – {tag}")
    plt.tight_layout()
    plt.savefig(OUT_DIR / f"{tag}_dotplot.png", dpi=300)
    plt.close()
    print(f"✓ done {tag}")


    # combined figure with heatmap (left) and dot‑plot (right)
    fig, (ax_hm, ax_dp) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 1.2]})
    
    # --- heatmap (ax_hm) ---
    im = ax_hm.imshow(mat, cmap="coolwarm", vmin=colors.min(), vmax=colors.max())
    ax_hm.set_xticks(np.arange(len(mat.columns)))
    ax_hm.set_xticklabels(mat.columns, rotation=90)
    ax_hm.set_yticks(np.arange(len(mat.index)))
    ax_hm.set_yticklabels(mat.index)
    for i in range(len(mat.index)):
        for j in range(len(mat.columns)):
            val = mat.iloc[i, j]
            if not np.isnan(val):
                ax_hm.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6)
    fig.colorbar(im, ax=ax_hm)
    ax_hm.set_title("Heatmap")

    # --- dot‑plot (ax_dp) ---
    scat = ax_dp.scatter(xpos, ypos, c=colors, s=sizes, cmap="jet", alpha=0.8)
    scat.set_clim(colors.min(), colors.max())
    cb = fig.colorbar(scat, ax=ax_dp, fraction=0.046, pad=0.04)
    cb.set_label("Spearman r", rotation=270, labelpad=10)
    for p in [0.05, 0.01, 0.001]:
        ax_dp.scatter([], [], s=20 * -np.log10(p + 1e-16), color="gray", alpha=0.8, label=f"p ≤ {p}")
    ax_dp.legend(title="p-value", loc="lower right", frameon=False)
    ax_dp.set_xticks(range(len(uniq_B)))
    ax_dp.set_xticklabels(uniq_B, rotation=90)
    ax_dp.set_yticks(range(len(uniq_A)))
    ax_dp.set_yticklabels(uniq_A)
    ax_dp.set_xlabel(f"stage – {ds2}")
    ax_dp.set_ylabel(f"stage – {ds1}")
    ax_dp.set_title("Dot‑plot")

    fig.suptitle(f"Stage‑level Spearman Correlation – {tag}")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(OUT_DIR / f"{tag}_combined.png", dpi=300)
    plt.close(fig)

    print(f"✓ done {tag} (saved combined figure)")


# ────── Run all three pairwise comparisons ──────────────────────────────────
for d1, d2 in combinations(DATASETS, 2):
    compute_pair(d1, d2)
print("All pairwise comparisons finished ✔")

