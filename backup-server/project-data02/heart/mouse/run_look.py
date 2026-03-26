#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_look.py
-----------
Read an AnnData (.h5ad), print a concise report about .obs, .var, and the whole object to a text file,
and randomly sample up to 20k cells to save as a subset .h5ad.

Usage:
  python run_look.py \
    --in ../data/GSE247719_PanSci_03_Heart_adata.h5ad \
    --out-info look_info.txt \
    --out-sub ../data/GSE247719_PanSci_03_Heart_adata.sub20k.h5ad \
    --n 20000 \
    --seed 42
"""
import argparse
import os
import sys
import gc
from datetime import datetime

try:
    import scanpy as sc
    import anndata as ad
    import numpy as np
    import pandas as pd
except Exception as e:
    print("Failed to import required packages. Make sure scanpy, anndata, numpy, pandas are installed.", file=sys.stderr)
    raise

def human_bytes(n: int) -> str:
    """Return human-readable bytes."""
    units = ["B","KB","MB","GB","TB"]
    size = float(n)
    for u in units:
        if size < 1024 or u == units[-1]:
            return f"{size:.2f} {u}"
        size /= 1024

def adata_memory_estimate(adata: "ad.AnnData") -> int:
    """Very rough memory footprint estimate (in bytes)."""
    total = 0
    # X
    X = adata.X
    try:
        total += X.data.nbytes + X.indptr.nbytes + X.indices.nbytes  # sparse
    except AttributeError:
        try:
            total += X.nbytes  # dense
        except Exception:
            pass
    # layers
    for k, M in (adata.layers.items() if adata.layers is not None else []):
        try:
            total += getattr(M, "nbytes", 0) or (M.data.nbytes + M.indptr.nbytes + M.indices.nbytes)
        except Exception:
            pass
    # obs/var
    for df in [adata.obs, adata.var]:
        try:
            total += df.memory_usage(deep=True).sum()
        except Exception:
            pass
    # obsm/varm/obsp/varp/uns (skip detailed recursion; just best-effort on arrays)
    def add_if_arraylike(obj):
        nonlocal total
        try:
            total += getattr(obj, "nbytes", 0)
        except Exception:
            try:
                total += obj.values.nbytes  # pandas
            except Exception:
                pass
    for mapping in [adata.obsm, adata.varm]:
        try:
            for v in mapping.values():
                add_if_arraylike(v)
        except Exception:
            pass
    for mapping in [adata.obsp, adata.varp]:
        try:
            for v in mapping.values():
                add_if_arraylike(v)
        except Exception:
            pass
    # uns (skip: often nested/python objects; can be very large/non-array)
    return int(total)

def write_report(adata: "ad.AnnData", out_txt: str):
    with open(out_txt, "w", encoding="utf-8") as f:
        def w(line=""):
            f.write(str(line) + "\n")

        w("# AnnData quick report")
        w(f"Generated: {datetime.now().isoformat(timespec='seconds')}")
        w("-" * 80)
        w("BASIC INFO")
        w("-" * 80)
        w(repr(adata))
        w(f"n_obs (cells): {adata.n_obs}")
        w(f"n_vars (genes/features): {adata.n_vars}")
        try:
            w(f"n_layers: {len(adata.layers) if adata.layers is not None else 0}")
            if adata.layers is not None:
                w(f"layers: {list(adata.layers.keys())}")
        except Exception:
            pass
        try:
            w(f"obsm keys: {list(adata.obsm.keys())}")
        except Exception:
            pass
        try:
            w(f"varm keys: {list(adata.varm.keys())}")
        except Exception:
            pass
        try:
            w(f"obsp keys: {list(adata.obsp.keys())}")
        except Exception:
            pass
        try:
            w(f"varp keys: {list(adata.varp.keys())}")
        except Exception:
            pass
        try:
            w(f"uns keys (top-level): {list(adata.uns.keys())[:30]}{' ...' if len(adata.uns.keys())>30 else ''}")
        except Exception:
            pass
        try:
            mem = adata_memory_estimate(adata)
            w(f"Estimated in-memory footprint: {human_bytes(mem)}")
        except Exception:
            pass

        w("")
        w("-" * 80)
        w(".obs OVERVIEW")
        w("-" * 80)
        try:
            w(f".obs shape: {adata.obs.shape}")
            w("Column names:")
            w(", ".join(map(str, adata.obs.columns.tolist())) or "(none)")
            w("")
            w("dtypes:")
            w(adata.obs.dtypes.to_string())
            w("")
            w("Head (first 10 rows):")
            w(adata.obs.head(10).to_string())
            w("")
            # value counts preview for up to 15 categorical columns
            cat_cols = [c for c, dt in adata.obs.dtypes.items() if str(dt) == "category"]
            for c in cat_cols[:15]:
                w(f"value_counts for obs['{c}'] (top 20):")
                w(adata.obs[c].value_counts(dropna=False).head(20).to_string())
                w("")
        except Exception as e:
            w(f"(Failed to summarize .obs: {e})")

        w("-" * 80)
        w(".var OVERVIEW")
        w("-" * 80)
        try:
            w(f".var shape: {adata.var.shape}")
            w("Column names:")
            w(", ".join(map(str, adata.var.columns.tolist())) or "(none)")
            w("")
            w("dtypes:")
            w(adata.var.dtypes.to_string())
            w("")
            w("Head (first 10 rows):")
            # .var can be very large; keep it small
            w(adata.var.head(10).to_string())
        except Exception as e:
            w(f"(Failed to summarize .var: {e})")

def main():
    parser = argparse.ArgumentParser(description="Inspect an AnnData and sample up to N cells.")
    parser.add_argument("--in", dest="in_path", type=str, required=False,
                        default="../data/GSE247719_PanSci_03_Heart_adata.h5ad",
                        help="Path to input .h5ad")
    parser.add_argument("--out-info", dest="out_info", type=str, required=False,
                        default="look_info.txt",
                        help="Path to write the info/report text file")
    parser.add_argument("--out-sub", dest="out_sub", type=str, required=False,
                        default="../data/GSE247719_PanSci_03_Heart_adata.sub20k.h5ad",
                        help="Path to write the subset .h5ad")
    parser.add_argument("--n", dest="n_cells", type=int, default=20000, help="Number of cells to sample (max)")
    parser.add_argument("--seed", dest="seed", type=int, default=42, help="Random seed for sampling")
    args = parser.parse_args()

    in_path = args.in_path
    out_info = args.out_info
    out_sub = args.out_sub
    n_cells = max(1, int(args.n_cells))
    seed = int(args.seed)

    if not os.path.exists(in_path):
        print(f"Input file not found: {in_path}", file=sys.stderr)
        sys.exit(1)

    # Read full object (backed='r' can be used for very large files, but we'll load into memory for easier ops)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Reading: {in_path}")
    adata = sc.read_h5ad(in_path, backed=None)
    print(f"[OK] AnnData loaded: {adata.n_obs} cells x {adata.n_vars} vars")

    # Write report
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Writing report to: {out_info}")
    write_report(adata, out_info)
    print("[OK] Report written.")

    # Sample cells (without replacement)
    n = min(n_cells, adata.n_obs)
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=n, replace=False)
    idx.sort()  # keep ascending index for better write locality

    # Subset (slice rows)
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Creating subset of {n} cells (seed={seed})")
    adata_sub = adata[idx, :].copy()

    # Save subset
    out_dir = os.path.dirname(out_sub)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    print(f"[{datetime.now().isoformat(timespec='seconds')}] Saving subset to: {out_sub}")
    adata_sub.write_h5ad(out_sub, compression="gzip")
    print("[OK] Subset saved.")

    # Clean up
    del adata_sub
    del adata
    gc.collect()
    print(f"[{datetime.now().isoformat(timespec='seconds')}] Done.")

if __name__ == "__main__":
    main()
