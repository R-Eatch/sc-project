#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hard-coded basal extraction + violin plots for 6 datasets.

Requirements
- Each .h5ad has:
  * obs['newcelltype']
  * layers['count']  (raw counts)

What it does
1) Read 6 adata from: ../{dataset}/1.subset/{dataset}_cleaned.h5ad
2) Subset basal cells by obs['newcelltype'] == 'basal' (case-insensitive)
3) Reset X to counts: X = layers['count'].copy()
4) Rename basal label to: {dataset}-Basal (stored back into obs['newcelltype'])
5) Concatenate all basal subsets
6) Normalize + log1p (for visualization)
7) Plot violin for genes: Krt5 Lcn2 Trp63 Acta2, grouped by {dataset}-Basal

Outputs (default)
- ./basal_violin_out/basal_all_6datasets.h5ad
- ./basal_violin_out/violin_basal_4genes_by_datasetBasal.png
- ./basal_violin_out/violin_{gene}_by_datasetBasal.png (per gene)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scanpy as sc
import anndata as ad

# -------------------------
# User settings (edit here)
# -------------------------
DATASETS = ["M-MG", "R-MG", "S-MG", "S-AG", "R-AG", "R-CG"]
PATH_TMPL = "../{dataset}/1.subset/{dataset}_cleaned.h5ad"
CELLTYPE_COL = "newcelltype"
COUNT_LAYER = "counts"
GENES = ["Krt5",'Acta2',"Lcn2", "Trp63"]
OUTDIR = Path("basal_violin_out")
TARGET_SUM = 1e4
DO_LOG1P = True
DPI = 300


def extract_basal(adata: ad.AnnData, dataset: str) -> ad.AnnData:
    if CELLTYPE_COL not in adata.obs:
        raise KeyError(
            f"{dataset}: obs['{CELLTYPE_COL}'] not found. Available obs columns: {list(adata.obs.columns)[:30]}"
        )
    if COUNT_LAYER not in adata.layers:
        raise KeyError(
            f"{dataset}: layers['{COUNT_LAYER}'] not found. Available layers: {list(adata.layers.keys())}."
        )

    mask = adata.obs[CELLTYPE_COL].astype(str).str.lower().eq("basal")
    basal = adata[mask].copy()
    if basal.n_obs == 0:
        raise ValueError(f"{dataset}: no basal cells found (obs['{CELLTYPE_COL}'] == 'basal').")

    # Reset X to raw counts
    basal.X = basal.layers[COUNT_LAYER].copy()

    # Rename label
    basal.obs[CELLTYPE_COL] = f"{dataset}-Basal"

    # Helpful extra column
    basal.obs["dataset"] = dataset

    return basal


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    basal_list: list[ad.AnnData] = []

    for ds in DATASETS:
        fp = Path(PATH_TMPL.format(dataset=ds))
        if not fp.exists():
            raise FileNotFoundError(f"Missing file for {ds}: {fp}")

        print(f"[read] {ds}: {fp}")
        a = sc.read_h5ad(fp)
        del a.obs['cellid']
        b = extract_basal(a, ds)
        basal_list.append(b)
        print(f"  basal cells: {b.n_obs:,} | genes: {b.n_vars:,}")

    # Concatenate basal-only objects
    # join='outer' makes it robust if var_names differ slightly across datasets.
    basal_all = ad.concat(
        basal_list,
        join="outer",
        label="dataset",
        keys=DATASETS,
        index_unique="-",
        merge="same",
        fill_value=0,
    )

    # Ensure numeric matrix
    if hasattr(basal_all.X, "dtype") and not np.issubdtype(basal_all.X.dtype, np.number):
        basal_all.X = basal_all.X.astype(np.float32)

    # Normalize for visualization
    sc.pp.normalize_total(basal_all, target_sum=TARGET_SUM)
    if DO_LOG1P:
        sc.pp.log1p(basal_all)

    # Save combined basal object
    out_h5ad = OUTDIR / "basal_all_6datasets.h5ad"
    basal_all.write(out_h5ad)
    print(f"[write] {out_h5ad}")

    # Filter genes that exist
    genes_present = [g for g in GENES if g in basal_all.var_names]
    genes_missing = [g for g in GENES if g not in basal_all.var_names]
    if len(genes_present) == 0:
        raise ValueError(
            f"None of the genes are present in var_names. Missing: {genes_missing}. "
            "Check var_names (symbol/case) or edit GENES."
        )
    if genes_missing:
        print(f"[warn] missing genes (skipped): {genes_missing}")

    # Plot
    import matplotlib.pyplot as plt

    sc.settings.set_figure_params(dpi=DPI)

    fig_multi = OUTDIR / "violin_basal_4genes_by_datasetBasal.png"
    sc.pl.violin(
        basal_all,
        keys=genes_present,
        groupby=CELLTYPE_COL,  # contains {dataset}-Basal
        rotation=45,
        stripplot=True,
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.tight_layout()
    plt.savefig(fig_multi, dpi=DPI, bbox_inches="tight")
    plt.close()
    print(f"[plot] {fig_multi}")

    for g in genes_present:
        out = OUTDIR / f"violin_{g}_by_datasetBasal.png"
        sc.pl.violin(
            basal_all,
            keys=g,
            groupby=CELLTYPE_COL,
            rotation=45,
            stripplot=True,
            jitter=0.4,
            show=False,
        )
        plt.tight_layout()
        plt.savefig(out, dpi=DPI, bbox_inches="tight")
        plt.close()
        print(f"[plot] {out}")

    print("Done.")
    fig=sc.pl.stacked_violin(basal_all,var_names=genes_present, groupby=CELLTYPE_COL, dendrogram=False,return_fig=True,show=False)
    #fig.subtitle(f"{dataset} stacked violin", fontsize=16)
    fig.savefig(f"{OUTDIR}/basal_violin_all.png",dpi=300,bbox_inches="tight")
    fig.savefig(f"{OUTDIR}/Basal_violin_all.pdf",dpi=300,bbox_inches="tight")
    figpath = OUTDIR / "Basal_all_dotplot.png"
    dp = sc.pl.dotplot(
        basal_all,
        var_names=genes_present,
        groupby=CELLTYPE_COL,
        use_raw=False,
        standard_scale=None,
        show=False,
        return_fig=True
        )
    dp.savefig(figpath, dpi=DPI, bbox_inches="tight")



if __name__ == "__main__":
    main()

