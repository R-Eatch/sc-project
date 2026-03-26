#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Global stage‑level correlation analysis for multiple scRNA‑seq datasets
────────────────────────────────────────────────────────────────────────
* Workflow:
  1. Load all datasets (M-MG, R-MG, S-MG).
  2. Concatenate into one object (outer join).
  3. Preprocess jointly (Normalize -> Log1p -> HVG).
  4. Compute Mean Expression for every (Dataset, Stage) combination.
  5. Calculate All-vs-All Spearman correlations.
  6. Plot global Heatmap and Dotplot.

* Output:
  • Global_stage_correlation.csv
  • Global_heatmap.png
  • Global_dotplot.png
"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pathlib
import sys
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore")

# ────── User‑configurable section ────────────────────────────────────────────
DATASETS      = ["M-MG", "R-MG", "S-MG"]       # Dataset IDs (Order matters for plotting)
CELL_COL      = "newcelltype"                  # Stage/Cell-type column in .obs
N_TOP_GENES   = 2000                           # HVG count (Updated to 2000)
DATA_TEMPLATE = "../../{ds}/1.subset/{ds}_cleaned.h5ad"  # Path template
OUT_DIR       = pathlib.Path("global_stage_correlation")
# ─────────────────────────────────────────────────────────────────────────────

OUT_DIR.mkdir(parents=True, exist_ok=True)
print(f"Datasets to merge: {DATASETS}")
print(f"Output directory : {OUT_DIR.absolute()}")

# ────── Helper Functions ─────────────────────────────────────────────────────

def load_dataset(ds_id: str) -> sc.AnnData:
    """Load a cleaned AnnData object."""
    path = pathlib.Path(DATA_TEMPLATE.format(ds=ds_id))
    if not path.exists():
        sys.exit(f"✗ Error: Cannot find file at {path}")
    print(f"  → Loading {ds_id}...")
    ad = sc.read_h5ad(path)
    
    # Clean up unnecessary columns to save memory during merge
    if "cellid" in ad.obs.columns:
        ad.obs.drop(columns=["cellid"], inplace=True)
    
    # Label the dataset source
    ad.obs['dataset_origin'] = ds_id
    
    # Ensure raw counts are in .X for the next steps if .layers['counts'] exists
    if "counts" in ad.layers:
        ad.X = ad.layers["counts"].copy()
    
    return ad

# ────── Main Execution Flow ──────────────────────────────────────────────────

def main():
    # 1. Load and Concatenate
    print("\n[1/5] Loading and merging datasets...")
    ad_list = []
    for ds in DATASETS:
        ad_list.append(load_dataset(ds))
    
    # Concatenate: join='outer' keeps all genes initially (fill=0)
    adata_full = sc.concat(ad_list, join='outer', label='batch', keys=DATASETS, fill_value=0)
    
    print(f"  → Combined object shape: {adata_full.shape}")

    # 2. Joint Preprocessing
    print("\n[2/5] Preprocessing combined data...")
    # Ensure metadata column is string
    adata_full.obs[CELL_COL] = adata_full.obs[CELL_COL].astype(str)
    
    # Create a unique group identifier: "Dataset | Stage"
    # This distinguishes "Stage1" in M-MG from "Stage1" in R-MG
    # This addresses the naming overlap issue by compounding the keys
    adata_full.obs['unique_group'] = (
        adata_full.obs['dataset_origin'].astype(str) + " | " + 
        adata_full.obs[CELL_COL].astype(str)
    )

    # Standard normalization flow
    sc.pp.normalize_total(adata_full)
    sc.pp.log1p(adata_full)
    
    # Select HVGs based on the joint variation
    print(f"  → Selecting top {N_TOP_GENES} HVGs...")
    sc.pp.highly_variable_genes(adata_full, n_top_genes=N_TOP_GENES, batch_key="dataset_origin")
    
    # Subset to HVGs
    adata_hvg = adata_full[:, adata_full.var.highly_variable].copy()
    print(f"  → Shape after HVG subset: {adata_hvg.shape}")

    # 3. Compute Mean Expression per Group (Pseudo-bulk)
    print("\n[3/5] Computing mean expression per group...")
    
    # Get sorted list of unique groups to ensure consistent order in matrix
    # We sort primarily by Dataset order defined in CONFIG, then by Stage name
    all_groups = []
    for ds in DATASETS:
        # Filter by dataset first to maintain dataset order in the plot axis
        stages = sorted(adata_hvg[adata_hvg.obs['dataset_origin'] == ds].obs[CELL_COL].unique())
        for stage in stages:
            # Reconstruct the key to match 'unique_group' format
            all_groups.append(f"{ds} | {stage}")
    
    mean_vectors = []
    valid_groups = []

    for group in all_groups:
        # Extract cells for this specific dataset-stage combination
        sub = adata_hvg[adata_hvg.obs['unique_group'] == group]
        
        if sub.n_obs < 3:
            print(f"  ⚠ Warning: Group '{group}' has <3 cells. Skipping.")
            continue
            
        # Compute mean across cells (axis 0)
        # Using numpy array conversion for safety
        if isinstance(sub.X, np.ndarray):
            mean_expr = sub.X.mean(axis=0)
        else:
            mean_expr = sub.X.mean(axis=0).A1  # .A1 flattens sparse matrix result
            
        mean_vectors.append(mean_expr)
        valid_groups.append(group)

    mean_matrix = np.array(mean_vectors) # Shape: (n_groups, n_genes)
    
    # 4. Compute Correlation Matrix
    print("\n[4/5] Calculating global correlation matrix...")
    n_groups = len(valid_groups)
    corr_mat = np.zeros((n_groups, n_groups))
    pval_mat = np.zeros((n_groups, n_groups))

    # Loop for symmetric calculation
    # Using simple loops is fast enough for ~30x30 groups
    for i in range(n_groups):
        for j in range(n_groups):
            vec_a = mean_matrix[i]
            vec_b = mean_matrix[j]
            
            rho, p = spearmanr(vec_a, vec_b, nan_policy='omit')
            corr_mat[i, j] = rho
            pval_mat[i, j] = p

    # Save CSV Data
    print("  → Saving CSV...")
    df_rows = []
    for i in range(n_groups):
        for j in range(n_groups):
            df_rows.append({
                "Group_A": valid_groups[i],
                "Group_B": valid_groups[j],
                "Correlation": corr_mat[i, j],
                "P_value": pval_mat[i, j]
            })
    pd.DataFrame(df_rows).to_csv(OUT_DIR / "Global_stage_correlation.csv", index=False)

    # 5. Visualization (One Chart Strategy)
    print("\n[5/5] Generating visualizations...")
    
    # Figure Size logic: scale with number of groups to avoid crowding
    fig_size = max(8, n_groups * 0.4) 
    
    # --- A. HEATMAP ---
    fig, ax = plt.subplots(figsize=(fig_size, fig_size * 0.9))
    im = ax.imshow(corr_mat, cmap="coolwarm", vmin=0.2, vmax=0.9) # Adjust vmin/vmax based on data range
    
    # Axis labels
    ax.set_xticks(np.arange(n_groups))
    ax.set_yticks(np.arange(n_groups))
    ax.set_xticklabels(valid_groups, rotation=90, fontsize=8)
    ax.set_yticklabels(valid_groups, fontsize=8)
    
    # Add values text
    for i in range(n_groups):
        for j in range(n_groups):
            val = corr_mat[i, j]
            # Only show text if correlation is decent to avoid visual clutter 
            color = "white" if abs(val) > 0.75 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=6)

    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman Correlation")
    plt.title(f"Global Stage Correlation\n({', '.join(DATASETS)})")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Global_heatmap.png", dpi=300)
    plt.savefig(OUT_DIR / "Global_heatmap.pdf", dpi=300)
    plt.close()

    # --- B. DOTPLOT (Matrix Style) ---
    # Convert matrix indices to list for scattering
    x_indices, y_indices = np.meshgrid(np.arange(n_groups), np.arange(n_groups))
    x_flat = x_indices.flatten()
    y_flat = y_indices.flatten()
    c_flat = corr_mat.flatten()
    p_flat = pval_mat.flatten()
    
    # Calculate sizes based on p-value (-log10)
    # Clip small p-values to avoid infinite size
    p_log = -np.log10(p_flat + 1e-20)
    sizes = 10 * p_log  # Scale factor
    
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))
    
    scat = ax.scatter(x_flat, y_flat, c=c_flat, s=sizes, cmap="jet", alpha=0.8, vmin=0.2, vmax=0.9)
    
    # Axes
    ax.set_xticks(np.arange(n_groups))
    ax.set_yticks(np.arange(n_groups))
    ax.set_xticklabels(valid_groups, rotation=90, fontsize=8)
    ax.set_yticklabels(valid_groups, fontsize=8)
    
    # Add grid to separate datasets visually if possible, or just standard grid
    ax.grid(True, linestyle='--', alpha=0.3)

    # Colorbar
    cbar = plt.colorbar(scat, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Spearman Correlation", rotation=270, labelpad=15)
    
    # P-value Legend
    # Create dummy points for legend
    legend_pvals = [0.05, 0.001, 1e-10]
    legend_handles = []
    for p in legend_pvals:
        s = 10 * -np.log10(p + 1e-20)
        legend_handles.append(ax.scatter([], [], s=s, color='gray', alpha=0.6, label=f'p < {p}'))
    
    ax.legend(handles=legend_handles, title="Significance", bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.title("Global Stage Correlation Dotplot")
    plt.tight_layout()
    plt.savefig(OUT_DIR / "Global_dotplot.png", dpi=300)
    plt.savefig(OUT_DIR / "Global_dotplot.pdf", dpi=300)
    plt.close()

    print(f"✓ Analysis Complete. Check output in: {OUT_DIR}")

if __name__ == "__main__":
    main()
