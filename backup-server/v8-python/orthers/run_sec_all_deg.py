#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Merge sec-like subtypes into LumSEC (per-dataset mapping) and run DEG per dataset.

User constraints (applied)
- Do NOT concatenate / merge adatas across datasets (each dataset processed independently)
- Do NOT reset X from layers['count']
- Do NOT run normalize_total / log1p (use the matrix as-is)

Datasets (fixed)
- M-MG, R-MG, S-MG, S-AG, R-AG, R-CG

Input paths (fixed)
- ../{dataset}/1.subset/{dataset}_cleaned.h5ad

Cell type column
- obs['subtype']  (edit SUBTYPE_COL if yours differs)

What this script does (per dataset)
1) Read adata
2) Create obs['subtype_merged']:
   - If obs['subtype'] is in SEC_MAP[dataset] => 'LumSEC'
   - Else keep original subtype
3) Run DEG with:
   sc.tl.rank_genes_groups(ea, groupby='subtype_merged', method='wilcoxon')
4) Export:
   - Top5 dotplot
   - DEG table via sc.get.rank_genes_groups_df(..., pval_cutoff=0.01, log2fc_min=0.25)

Outputs
- ./sec_merge_deg_out/{dataset}/{dataset}_ranked_genes_for_subtype_merged.csv
- ./sec_merge_deg_out/{dataset}/{dataset}_dotplot_subtype_merged_top5.png
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc

# -------------------------
# Fixed configuration
# -------------------------
DATASETS = ["M-MG", "R-MG", "S-MG", "S-AG", "R-AG", "R-CG"]
PATH_TMPL = "../{dataset}/1.subset/{dataset}_cleaned.h5ad"
SUBTYPE_COL = "newcelltype"          # original subtype column (edit if needed)
MERGED_COL = "subtype_merged"    # new merged column
OUTDIR = Path("sec_merge_deg_out")

# DEG params (match your reference)
METHOD = "wilcoxon"
PVAL_CUTOFF = 0.01
LOG2FC_MIN = 0.25
N_GENES_PLOT = 5
MIN_LOGFOLCHANGE_PLOT = 0.25
DPI = 300

# Per-dataset mapping: which original subtypes should be merged into LumSEC
SEC_MAP: dict[str, set[str]] = {
    "M-MG": {"LumSEC-Lac", "LumSEC-Lip"},
    "R-MG": {"LumSEC-Lac", "LumSEC-Lip", "LumSEC-Mgp", "LumSEC-Vacm1"},
    "S-MG": {"LumSEC-Lac", "LumSEC-Lip"},
    "R-AG": {"LumSEC-AG-t1", "LumSEC-AG-t2"},
    "R-CG": {"LumSEC-Lip-CG", "Lum-Tm4sf4", "Lum-Stat4"},
    "S-AG": {"LumSEC-AG-Pip"},
}


def merge_sec_to_lumsec(adata, dataset: str):
    if SUBTYPE_COL not in adata.obs:
        raise KeyError(
            f"{dataset}: obs['{SUBTYPE_COL}'] not found. Available obs columns: {list(adata.obs.columns)[:40]}"
        )

    sec_set = SEC_MAP.get(dataset)
    if sec_set is None:
        raise KeyError(f"No SEC_MAP entry for dataset '{dataset}'.")

    sub = adata.obs[SUBTYPE_COL].astype(str)
    adata.obs[MERGED_COL] = np.where(sub.isin(sec_set), "LumSEC", sub)

    n_merged = int((adata.obs[MERGED_COL].values == "LumSEC").sum())
    print(f"  merged to LumSEC: {n_merged:,} cells")

    return adata


def do_DEG_merged(ea, dataset: str, outdir: Path):
    """Run DEG using merged subtype column and export results + top5 dotplot.

    Important: uses ea.X as-is (no layer reset, no normalization/log1p).
    """

    # Ensure groupby col is categorical for stable ordering
    if ea.obs[MERGED_COL].dtype.name != "category":
        ea.obs[MERGED_COL] = ea.obs[MERGED_COL].astype("category")

    # DEG
    sc.tl.rank_genes_groups(ea, groupby=MERGED_COL, method=METHOD)

    # Plot
    import matplotlib.pyplot as plt

    figpath = outdir / f"{dataset}_dotplot_subtype_merged_top{N_GENES_PLOT}.png"
    sc.pl.rank_genes_groups_dotplot(
        ea,
        groupby=MERGED_COL,
        n_genes=N_GENES_PLOT,
        min_logfoldchange=MIN_LOGFOLCHANGE_PLOT,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(figpath, dpi=DPI, bbox_inches="tight")
    plt.close()

    # Export DEG table
    deg_df = sc.get.rank_genes_groups_df(
        ea,
        group=None,
        pval_cutoff=PVAL_CUTOFF,
        log2fc_min=LOG2FC_MIN,
    )
    out_csv = outdir / f"{dataset}_ranked_genes_for_subtype_merged.csv"
    deg_df.to_csv(out_csv, index=False)

    print(f"  [write] {out_csv}")
    print(f"  [plot ] {figpath}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for ds in DATASETS:
        fp = Path(PATH_TMPL.format(dataset=ds))
        if not fp.exists():
            raise FileNotFoundError(f"Missing file for {ds}: {fp}")

        ds_out = OUTDIR / ds
        ds_out.mkdir(parents=True, exist_ok=True)

        print(f"[read] {ds}: {fp}")
        adata = sc.read_h5ad(fp)

        # Merge sec-like subtypes into LumSEC (only changes obs column)
        adata = merge_sec_to_lumsec(adata, ds)

        # Run DEG using ea.X as-is
        do_DEG_merged(adata, ds, ds_out)

    print("Done.")


if __name__ == "__main__":
    main()

