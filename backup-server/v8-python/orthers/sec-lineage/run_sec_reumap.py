#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-subset single-cell datasets by `newcelltype`, re-run normalization + UMAP,
then export rasterized (gridized) UMAP PDFs; additionally plot per-dataset genelist
heatmaps across newcelltype using adata.layers['normalized'].

User rules
- Exclude (case-sensitive, if present): Basal, LumHR, Epi-Fibro, Lum-Basal, Lum-immune, Lum-Immune
- Re-normalize after subsetting, then HVG -> PCA(PC=20) -> neighbors -> leiden(res=0.3) -> UMAP
- Plot UMAP using sc.pl.umap(color='newcelltype') and rasterize the scatter layer, save as PDF
- Title of UMAP figure: dataset only
- Save processed (filtered + re-UMAP) AnnData as .h5ad
- Heatmap per dataset: genes (rows) x newcelltype (cols) using adata.layers['normalized']
  - Aggregate: mean expression per newcelltype
  - Transform for display: per-gene z-score across newcelltype (configurable)
  - Export as PDF, red/blue colormap

Outputs
- result/reumap_newcelltype/<dataset>/<dataset>_umap_newcelltype_grid.pdf
- result/reumap_newcelltype/<dataset>/<dataset>_filtered_reumap.h5ad
- result/reumap_newcelltype/<dataset>/<dataset>_heatmap_genelist_newcelltype.pdf

Requirements
- scanpy, anndata, numpy, pandas, matplotlib, harmonypy
"""

from __future__ import annotations

from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import harmonypy as hm

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# =========================
# Global configuration
# =========================
DATASETS = ["M-MG", "R-MG", "S-MG", "S-AG", "R-AG", "R-CG"]
PATH_TMPL = "../../{dataset}/1.subset/{dataset}_cleaned.h5ad"

# Per-dataset path override (user request)
OVERRIDE_PATHS = {
    "S-MG": "/data01/sunxuebo/project/scrnaseq/v8-python/orthers/S-MG-holo.h5ad",
}

NEWCELLTYPE_COL = "newcelltype"

# Case-sensitive exclusions
EXCLUDE_CELLTYPES = [
    "Basal",
    "LumHR",
    "Epi-Fibro",
    "Lum-Basal",
    "Lum-immune",
    "Lum-Immune",
]

# Per-dataset secondary exclusions (applied AFTER EXCLUDE_CELLTYPES)
SECONDARY_EXCLUDE_BY_DATASET = {
    # After the first round of filtering, drop extra celltypes for specific datasets
    "S-MG": ["Epi-Lgals7"],
    "R-AG": ["Epi-Krt7", "Epi-Lgals7", "Epi-Pro"],
}

TOP_GENES = 2000
PC = 10
RESOLUTION = 0.3

RUN_HARMONY = True
USE_SCVI = False

RANDOM_STATE = 42
N_NEIGHBORS = 15

OUTROOT = Path("result/reumap_newcelltype")
DPI = 300

# Heatmap settings
GENELIST = [
    # TODO: put your genes here
    # "Lcn2", "Trp63", "Krt14",
]
HEATMAP_LAYER = "normalized"  # user request
HEATMAP_ZSCORE_PER_GENE = False
HEATMAP_FIGSIZE = (11, 0.6)  # wider and less flat; height is scaled by genes  # (width, height_per_gene) -> auto-height
HEATMAP_CMAP = "bwr"  # red/blue
HEATMAP_VLIM = 2.5  # symmetric for z-score display

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120, facecolor="white")


# =========================
# User-provided functions (kept as-is, minimal changes)
# =========================
def run_preprocess(adata, top_genes: int):
    # ensure counts layer exists
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    adata.layers["normalized"] = adata.X.copy()
    print("Finish normalized")
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
    sc.pl.highly_variable_genes(adata, show=False)
    print("Finish Varible genes")


def run_harmony(adata, vars_use: str = "sample"):
    print("running harmony")
    pca_result = adata.obsm["X_pca"]
    ho = hm.run_harmony(pca_result, adata.obs, vars_use, random_state=42)
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T
    print("finished harmony")
    return adata


def run_reduceDimension(ea, use_scvi: bool, runharmony: bool, PCs: int, res: float):
    sc.pp.scale(ea, max_value=10)
    sc.tl.pca(ea, mask_var="highly_variable")
    sc.pl.pca_variance_ratio(ea, log=False, show=False)

    if runharmony:
        run_harmony(ea)
        sc.pp.neighbors(
            ea,
            use_rep="X_pca_harmony",
            n_pcs=PCs,
            random_state=RANDOM_STATE,
            n_neighbors=N_NEIGHBORS,
        )
        print("Finish neighbors")
    else:
        sc.pp.neighbors(
            ea,
            n_pcs=PCs,
            random_state=RANDOM_STATE,
            n_neighbors=N_NEIGHBORS,
        )
        print("Finish neighbors")

    sc.tl.leiden(ea, resolution=res, random_state=RANDOM_STATE)
    print("Finish clustering")
    sc.tl.umap(ea, random_state=RANDOM_STATE)
    print("Finish UMAP")


# =========================
# Core helpers
# =========================
def subset_by_exclusion(adata, group_col: str, exclude_labels: List[str]):
    """Case-sensitive exclusion by exact label match in obs[group_col]."""
    if group_col not in adata.obs.columns:
        raise KeyError(
            f"obs column '{group_col}' not found. Available: {list(adata.obs.columns)[:30]} ..."
        )

    s = adata.obs[group_col]

    present = set(s.dropna().unique().tolist())
    actually_excluded = sorted(list(present.intersection(set(exclude_labels))))

    keep_mask = ~s.isin(exclude_labels)
    sub = adata[keep_mask].copy()
    return sub, actually_excluded


def rasterize_scanpy_scatter(fig: matplotlib.figure.Figure):
    """Rasterize PathCollection(s) in a Scanpy plot, keeping text/legend vector."""
    for ax in fig.axes:
        for col in ax.collections:
            try:
                col.set_rasterized(True)
            except Exception:
                pass


def plot_umap_gridpdf_scanpy(adata, group_col: str, out_pdf: Path, title: str):
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    sc.pl.umap(
        adata,
        color=group_col,
        title=title,
        show=False,
    )

    fig = plt.gcf()
    rasterize_scanpy_scatter(fig)

    fig.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def _to_dense_1d(x):
    """Convert 1D slice (n_cells,) to dense numpy array."""
    try:
        import scipy.sparse as sp

        if sp.issparse(x):
            return np.asarray(x.todense()).ravel()
    except Exception:
        pass
    return np.asarray(x).ravel()


def plot_genelist_heatmap_by_newcelltype(
    adata,
    genelist: List[str],
    group_col: str,
    layer: str,
    out_pdf: Path,
    zscore_per_gene: bool = True,
):
    """Heatmap of mean expression per newcelltype using a specified layer.

    Matrix construction:
    - For each gene in genelist (rows) and each newcelltype (cols), compute mean over cells
      using adata.layers[layer].

    Display transform:
    - Optionally z-score per gene across columns, so each row highlights relative enrichment
      across celltypes (not absolute expression differences between genes).
    """
    out_pdf.parent.mkdir(parents=True, exist_ok=True)

    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in adata.layers. Available: {list(adata.layers.keys())}")

    if group_col not in adata.obs.columns:
        raise KeyError(f"obs column '{group_col}' not found")

    # keep genes present
    present_genes = [g for g in genelist if g in adata.var_names]
    missing_genes = [g for g in genelist if g not in adata.var_names]
    if len(present_genes) == 0:
        print(f"[HEATMAP][SKIP] No genes found in var_names. Missing: {missing_genes[:20]}{'...' if len(missing_genes)>20 else ''}")
        return

    # order celltypes (stable by appearance)
    s = adata.obs[group_col]
    celltypes = pd.unique(s.dropna())
    celltypes = [ct for ct in celltypes if ct != ""]

    # index genes
    gene_idx = {g: i for i, g in enumerate(adata.var_names)}
    g_indices = np.array([gene_idx[g] for g in present_genes], dtype=int)

    X = adata.layers[layer]

    # build mean matrix: genes x celltypes
    mat = np.zeros((len(present_genes), len(celltypes)), dtype=float)

    for j, ct in enumerate(celltypes):
        mask = (s == ct).values
        if mask.sum() == 0:
            continue
        # slice: cells x genes
        X_sub = X[mask][:, g_indices]
        # mean over cells
        try:
            import scipy.sparse as sp

            if sp.issparse(X_sub):
                mean_vec = np.asarray(X_sub.mean(axis=0)).ravel()
            else:
                mean_vec = X_sub.mean(axis=0)
        except Exception:
            mean_vec = np.asarray(X_sub).mean(axis=0)
        mat[:, j] = np.asarray(mean_vec).ravel()

    disp = mat.copy()

    # z-score per gene for display (row-wise)
    if zscore_per_gene:
        mu = disp.mean(axis=1, keepdims=True)
        sd = disp.std(axis=1, keepdims=True)
        sd[sd == 0] = 1.0
        disp = (disp - mu) / sd

        # Plot
    height = max(4.0, HEATMAP_FIGSIZE[1] * len(present_genes))
    fig, ax = plt.subplots(figsize=(HEATMAP_FIGSIZE[0], height))

    if zscore_per_gene:
        vmin, vmax = -HEATMAP_VLIM, HEATMAP_VLIM
    else:
        vmin, vmax = np.percentile(disp, 1), np.percentile(disp, 99)

    im = ax.imshow(disp, aspect="auto", interpolation="nearest", cmap=HEATMAP_CMAP, vmin=vmin, vmax=vmax)

    ax.set_yticks(np.arange(len(present_genes)))
    ax.set_yticklabels(present_genes, fontsize=8)

    ax.set_xticks(np.arange(len(celltypes)))
    ax.set_xticklabels(celltypes, rotation=45, ha="right", fontsize=8)

    ax.set_xlabel(group_col)
    ax.set_ylabel("Genes")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Z-score" if zscore_per_gene else "Mean expression")

    ax.set_title(f"Heatmap ({layer})")

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


# =========================
# Main
# =========================
def main():
    OUTROOT.mkdir(parents=True, exist_ok=True)

    for ds in DATASETS:
        in_path_str = OVERRIDE_PATHS.get(ds, PATH_TMPL.format(dataset=ds))
        in_path = Path(in_path_str)

        if not in_path.exists():
            print(f"[SKIP] Missing file: {in_path}")
            continue

        print("" + "=" * 80)
        print(f"[LOAD] {ds}: {in_path}")
        adata = sc.read_h5ad(in_path)

        if NEWCELLTYPE_COL not in adata.obs.columns:
            print(
                f"[SKIP] '{NEWCELLTYPE_COL}' not found in {ds}. "
                f"Available obs cols: {list(adata.obs.columns)[:30]} ..."
            )
            continue

        # 1) subset by exclusion list
        before_n = adata.n_obs
        adata, excluded = subset_by_exclusion(adata, NEWCELLTYPE_COL, EXCLUDE_CELLTYPES)
        after_n = adata.n_obs
        print(f"[SUBSET] {ds}: {before_n} -> {after_n} cells")
        print(
            f"[EXCLUDED] {ds}: {excluded if excluded else 'None of the excluded labels were present (case-sensitive)'}"
        )

        # 1b) optional per-dataset secondary exclusion (after the first filter)
        secondary_exclude = SECONDARY_EXCLUDE_BY_DATASET.get(ds, [])
        if len(secondary_exclude) > 0:
            before_n2 = adata.n_obs
            adata, excluded2 = subset_by_exclusion(adata, NEWCELLTYPE_COL, secondary_exclude)
            after_n2 = adata.n_obs
            print(f"[SUBSET2] {ds}: {before_n2} -> {after_n2} cells")
            print(
                f"[EXCLUDED2] {ds}: {excluded2 if excluded2 else 'None of the secondary excluded labels were present (case-sensitive)'}"
            )

        # 2) re-normalize + HVG
        run_preprocess(adata, TOP_GENES)

        # 3) PCA / neighbors / leiden / umap
        #run_reduceDimension(adata, USE_SCVI, RUN_HARMONY, PCs=PC, res=RESOLUTION)

        # 4) plot rasterized UMAP PDF (via sc.pl.umap), title is dataset only
        out_pdf = OUTROOT / ds / f"{ds}_umap_{NEWCELLTYPE_COL}_grid.pdf"
        plot_umap_gridpdf_scanpy(adata, NEWCELLTYPE_COL, out_pdf, title=ds)
        print(f"[SAVE] {out_pdf}")

        # 5) save processed h5ad
        out_h5ad = OUTROOT / ds / f"{ds}_filtered_reumap.h5ad"
        out_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata.write(out_h5ad)
        print(f"[SAVE] {out_h5ad}")

        # 6) heatmap per dataset (if GENELIST is provided)
        if len(GENELIST) > 0:
            out_hm = OUTROOT / ds / f"{ds}_heatmap_genelist_{NEWCELLTYPE_COL}.pdf"
            plot_genelist_heatmap_by_newcelltype(
                adata,
                genelist=GENELIST,
                group_col=NEWCELLTYPE_COL,
                layer=HEATMAP_LAYER,
                out_pdf=out_hm,
                zscore_per_gene=HEATMAP_ZSCORE_PER_GENE,
            )
            print(f"[SAVE] {out_hm}")
        else:
            print(f"[HEATMAP][SKIP] GENELIST is empty; no heatmap generated for {ds}")

    print("Done.")


if __name__ == "__main__":
    main()

