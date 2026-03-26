#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_subset_by_obs.py
--------------------
Subset an AnnData by obs columns:
  - Age_group in {"03_months", "23_months"}
  - Main_cell_type in {"Fibroblasts", "Ventricular cardiomyocytes", "Endocardial endothelial cells"}

Reads:  ../data/GSE247719_PanSci_03_Heart_adata.h5ad
Writes: ../data/GSE247719_PanSci_03_Heart_adata.sub_Age03_23_Fibro_VentCM_EndoEndo.h5ad

Usage:
  python run_subset_by_obs.py
"""

import sys
import os
from datetime import datetime

try:
    import scanpy as sc
    import anndata as ad
    import pandas as pd
    import numpy as np
except Exception as e:
    print("Failed to import required packages. Make sure scanpy, anndata, pandas, numpy are installed.", file=sys.stderr)
    raise

# ---- Config (hard-coded as requested) ----
IN_PATH  = "../data/GSE247719_PanSci_03_Heart_adata.h5ad"
#OUT_PATH = "../data/GSE247719_PanSci_03_Heart_adata.sub_Age03_23_Fibro_VentCM_EndoEndo.h5ad"
OUT_PATH = "../data/GSE247719_PanSci_03_Heart_adata.sub_Age03_23_Main_cells.h5ad"
AGE_COL = "Age_group"
AGE_KEEP = ["03_months", "23_months"]

TYPE_COL = "Main_cell_type"
#TYPE_KEEP = ["Fibroblasts", "Ventricular cardiomyocytes", "Endocardial endothelial cells"]
TYPE_KEEP = [
    "Vascular endothelial cells",
    "Fibroblasts",
    "Ventricular cardiomyocytes",
    "Myeloid cells",
    "Mural cells",
    "Endocardial endothelial cells",
    "Lymphoid cells_T cells",
    "Lymphatic endothelial cells",
    "Atrial cardiomyocytes",
    "Lymphoid cells_B cells",
    "Mesothelial cells",
    "Neural cells",
    "Adipocytes"
]
def main():
    if not os.path.exists(IN_PATH):
        print(f"[ERROR] Input not found: {IN_PATH}", file=sys.stderr)
        sys.exit(1)

    print(f"[{datetime.now().isoformat(timespec='seconds')}] Reading: {IN_PATH}")
    adata = sc.read_h5ad(IN_PATH, backed=None)

    adata.X = adata.layers['raw_counts'].copy()
    print(f"[OK] Loaded AnnData: {adata.shape[0]} cells x {adata.shape[1]} vars")

    # ---- sanity checks ----
    for col in (AGE_COL, TYPE_COL):
        if col not in adata.obs.columns:
            print(f"[ERROR] Required obs column '{col}' not in adata.obs. Available: {list(adata.obs.columns)}", file=sys.stderr)
            sys.exit(2)

    # ---- build masks ----
    age_mask = adata.obs[AGE_COL].isin(AGE_KEEP)
    type_mask = adata.obs[TYPE_COL].isin(TYPE_KEEP)
    mask = age_mask & type_mask

    n_match = int(mask.sum())
    print(f"[INFO] Cells passing filters: {n_match} / {adata.n_obs}")
    if n_match == 0:
        print("[ERROR] No cells matched the given filters. Please check values in obs.", file=sys.stderr)
        # print unique values to help debugging
        print(f"Unique {AGE_COL}: {adata.obs[AGE_COL].astype(str).unique()[:50]}", file=sys.stderr)
        print(f"Unique {TYPE_COL}: {adata.obs[TYPE_COL].astype(str).unique()[:50]}", file=sys.stderr)
        sys.exit(3)

    # ---- subset ----
    adata_sub = adata[mask, :].copy()

    # drop unused categories (if any) for tidier object
    for df in (adata_sub.obs, adata_sub.var):
        for c in df.select_dtypes(include=["category"]).columns:
            df[c] = df[c].cat.remove_unused_categories()

    # ---- save ----
    out_dir = os.path.dirname(OUT_PATH)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"[{datetime.now().isoformat(timespec='seconds')}] Saving subset to: {OUT_PATH}")
    adata_sub.write_h5ad(OUT_PATH, compression="gzip")
    print("[OK] Done.")

if __name__ == "__main__":
    main()
