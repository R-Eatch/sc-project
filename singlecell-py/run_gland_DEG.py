#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
DEG pipeline (Scanpy):
- Rename glands:  SG -> EG
- Merge glands:   CG + AG -> SG  (after the rename above)
- Per-gland DEG with groupby=anno (Wilcoxon)
- Output: one dotplot PNG per gland, one CSV per gland (only pval/log2fc filtered results; NO top-N clipping)
- Robust to missing var['gene_name']

Assumptions:
- Current adata.X is already log1p-normalized (no normalization/log1p will be done here).

Required obs fields:
- adata.obs['gland']  (original may contain: MG, AG, SG, CG)
- adata.obs['anno']
"""

from pathlib import Path
import scanpy as sc
import pandas as pd

# ---------------- Config ----------------
INPUT_H5AD    = "../../h5ad2seurat/processed/v9_all_counts.h5ad"  # ← updated path
OUTDIR        = "deg_out"           # output dir
DATASET_TAG   = "MyDataset"         # for filenames
N_GENES_PLOT  = 5                    # genes per group in dotplot
PVALUE_CUTOFF = 0.01
LOG2FC_MIN    = 0.25

# ---------------- Utils ----------------
def ensure_dir(p: str | Path):
    Path(p).mkdir(parents=True, exist_ok=True)


def ensure_gene_name_var(adata):
    if 'gene_name' not in adata.var.columns:
        adata.var['gene_name'] = adata.var_names.astype(str)
        print("[FIX] var['gene_name'] 不存在，已用 var_names 填充。")


def rename_and_merge_glands(adata):
    """Apply:
      1) SG -> EG
      2) AG -> SG, CG -> SG  (after step 1)
    """
    if 'gland' not in adata.obs.columns:
        raise ValueError("缺少 adata.obs['gland'] 列。")

    # Step 1: SG -> EG
    adata.obs['gland'] = adata.obs['gland'].replace({'SG': 'EG'})
    # Step 2: merge AG/CG into SG
    adata.obs['gland'] = adata.obs['gland'].replace({'AG': 'SG', 'CG': 'SG'})

    print("[OK] 已重命名/合并 gland：SG→EG；AG+CG→SG。")
    print("[INFO] 重命名后 gland 计数：")
    print(adata.obs['gland'].value_counts())


def export_deg_per_gland(adata, outdir: str | Path, dataset_tag: str):
    ensure_dir(outdir)
    glands = sorted(adata.obs['gland'].astype(str).unique().tolist())
    print(f"[INFO] 检测到 gland：{glands}")

    for g in glands:
        print(f"[INFO] 处理 gland = {g}")
        ad = adata[adata.obs['gland'].astype(str) == g].copy()

        if 'anno' not in ad.obs.columns:
            raise ValueError("缺少 ad.obs['anno']。")
        if ad.obs['anno'].nunique() < 2:
            print(f"[WARN] gland={g} 的 anno 组别不足 2 个，跳过。")
            continue

        # Ensure var['gene_name'] exists in the subset too
        #ensure_gene_name_var(ad)

        # DEG (Wilcoxon)
        sc.tl.rank_genes_groups(
            ad,
            groupby='anno',
            method='wilcoxon',
            use_raw=False
        )

        # Dotplot — only pass gene_symbols if the column exists
        dotplot_kwargs = dict(
            adata=ad,
            groupby='anno',
            n_genes=N_GENES_PLOT,
            min_logfoldchange=LOG2FC_MIN,
            show=False,
            save=f"_{dataset_tag}_dotplot_{g}.png",
        )
        if 'gene_name' in ad.var.columns:
            dotplot_kwargs['gene_symbols'] = 'gene_name'
        sc.pl.rank_genes_groups_dotplot(**dotplot_kwargs)

        # DEG table — pass gene_symbols only if present
        get_df_kwargs = dict(
            adata=ad,
            group=None,
            pval_cutoff=PVALUE_CUTOFF,
            log2fc_min=LOG2FC_MIN,
        )
        if 'gene_name' in ad.var.columns:
            get_df_kwargs['gene_symbols'] = 'gene_name'
        deg_df = sc.get.rank_genes_groups_df(**get_df_kwargs)

        out_csv = Path(outdir) / f"{dataset_tag}_DEG_{g}.csv"
        deg_df.to_csv(out_csv, index=False, encoding='utf-8')
        print(f"[OK] 保存：{out_csv}")


# ---------------- Main ----------------
if __name__ == "__main__":
    sc.settings.verbosity = 2
    sc.settings.figdir = "figures"
    ensure_dir(sc.settings.figdir)
    ensure_dir(OUTDIR)

    print(f"[INFO] 读取：{INPUT_H5AD}")
    adata = sc.read_h5ad(INPUT_H5AD)

    # Ensure var['gene_name'] exists
    #ensure_gene_name_var(adata)

    # Rename & merge glands as requested
    rename_and_merge_glands(adata)

    # Export per-gland DEG
    export_deg_per_gland(adata, OUTDIR, DATASET_TAG)

    print("\n[DONE] 逐 gland DEG 导出完成。")

