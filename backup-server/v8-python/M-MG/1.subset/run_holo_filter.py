#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Scanpy two-dataset merge & embedding pipeline
--------------------------------------------
在 import 后可通过修改全局变量进行配置：
- 先对 adata1：按指定 leiden 列与集群列表进行子集化，并写入新列 newcelltype（可由映射或固定标签赋值）。
- 读取 adata2（其应已有 newcelltype 列）。
- 合并 adata1_subset 与 adata2；将 X 替换为 layers['normalized']（若缺失则自动计算并创建）。
- 标准流程：HVG → Scale → PCA → (可选 Harmony) → Neighbors → Leiden → UMAP。
- 输出：UMAP（按 ['newcelltype','stage','sample'] 着色）、FeaturePlot、cell barplot（stage/sample），并保存合并后的 h5ad。

依赖：scanpy, anndata, numpy, pandas, matplotlib, harmonypy（可选）。
"""

import os
os.environ["PYTHONHASHSEED"] = "0"

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 可选：Harmony（harmonypy）。若未安装，可将 RUN_HARMONY=False。
try:
    import harmonypy as hm
    _HAS_HARMONY = True
except Exception:
    _HAS_HARMONY = False

# =========================
# Global Variables（可修改）
# =========================
# 路径请按需修改
dataset = "M-MG"
ADATA1_PATH = f"{dataset}_for_DEG.h5ad"    # 需要被子集化的对象
ADATA2_PATH = f"{dataset}_cleaned.h5ad"    # 已包含 newcelltype 的对象
OUTPUT_DIR  = "./"
MERGED_NAME = f"{dataset}-holo.h5ad"    # 合并后保存的文件名
PREFIX      = f"{dataset}-holo"         # 输出图像前缀

# 元数据列名（如与实际不符请改）
LEIDEN_COL  = "leiden"
STAGE_COL   = "stage"
SAMPLE_COL  = "sample"
NEWCELL_COL = "newcelltype"
# 要从 .obs 中删除的列（可按需增减）
OBS_COLS_TO_DROP = ["cellid"]

# 需要保留的 adata1 集群（支持 int 或 str；内部统一转为 str 比较）
CLUSTERS_TO_KEEP = [5]

# 为 adata1 写入 newcelltype 的方式（2 选 1，优先使用映射）
# 方式 A：按集群映射到细胞类型名
# 形式 1：{"Sebocytes": [2,5], "Duct": [4], ...}
# 形式 2：{"2": "Sebocytes", "5": "Sebocytes", "4": "Duct"}
NEWCELLTYPE_MAPPING = {
    "Sebocytes": [5],
}
# 方式 B：为 adata1 的所有保留细胞赋一个固定标签（当上面的映射为空时生效）
FIXED_NEWCELLTYPE_LABEL = "Adata1_Selected"

# 处理参数
N_TOP_GENES  = 2000
N_PCS        = 20
N_NEIGHBORS  = 15
LEIDEN_RES   = 0.3
RUN_HARMONY  = True              # 若未安装 harmonypy 会被强制置为 False
HARMONY_VARS = [SAMPLE_COL]      # Harmony 按哪些 batch 变量校正
RANDOM_SEED  = 42

# FeaturePlot 基因
FEATURE_GENES = [
    "Plin2", "Cidea", "Fabp3", "Mgst1", "Acsl1", "Pdk4", "Xdh", "Tkt", "Echdc1",
    "Adipoq", "Awat2", "Epcam", "Esr1", "Lalba", "Krt5"
]

# 其他绘图参数
FIG_DPI  = 200
POINT_SZ = 3

# =========================
# Utility Functions
# =========================

def _ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)
    sc.settings.figdir = path


def _to_str_list(vals):
    return [str(v) for v in vals]

def drop_obs_columns(adata, cols):
    """Drop specified columns from adata.obs if present."""
    if not cols:
        return adata
    keep = [c for c in cols if c in adata.obs.columns]
    if keep:
        adata.obs.drop(columns=keep, inplace=True)
        print(f"[OBS] Dropped columns from obs: {keep}")
    else:
        print("[OBS] No target columns found in obs.")
    return adata

def ensure_normalized_layer(adata: ad.AnnData, layer_name: str = "normalized") -> ad.AnnData:
    """若不存在指定 normalized 层，则以 total-count normalize + log1p 计算并创建。
    保留原 X；新建/覆盖 adata.layers[layer_name]。
    """
    if layer_name in (adata.layers.keys() if adata.layers is not None else []):
        return adata
    # 以副本在内存中计算 normalized 表达（不改动 adata.X）
    tmp = adata.copy()
    sc.pp.normalize_total(tmp, target_sum=1e4)
    sc.pp.log1p(tmp)
    adata.layers[layer_name] = tmp.X.copy()
    return adata


def set_X_from_layer(adata: ad.AnnData, layer_name: str = "normalized") -> ad.AnnData:
    """将 adata.X 替换为指定层。若该层缺失，则自动创建后再替换。"""
    adata = ensure_normalized_layer(adata, layer_name)
    adata.X = adata.layers[layer_name].copy()
    return adata


def subset_and_set_newcelltype(
    adata: ad.AnnData,
    leiden_col: str,
    clusters_to_keep,
    newcelltype_mapping: dict | None = None,
    fixed_label: str | None = None,
    newcell_col: str = NEWCELL_COL,
) -> ad.AnnData:
    """对 adata 进行子集化，并写入 newcelltype 列。
    - 支持两种映射格式：
      1) {celltype: [clusters, ...]}
      2) {cluster: celltype}
    - 若 mapping 为空，则为所有保留细胞赋 fixed_label。
    """
    if leiden_col not in adata.obs:
        raise ValueError(f"'{leiden_col}' not in adata.obs; got columns: {list(adata.obs.columns)}")

    # 统一字符串比较
    keep_set = set(_to_str_list(clusters_to_keep))
    # 如果 leiden 是 category，先转为 str
    leiden_as_str = adata.obs[leiden_col].astype(str)
    subset_mask = leiden_as_str.isin(keep_set)
    adata = adata[subset_mask].copy()

    # 赋值 newcelltype
    if newcelltype_mapping:
        # 兼容两种映射格式
        mapping_cluster_to_type: dict[str, str] = {}
        for k, v in newcelltype_mapping.items():
            if isinstance(v, (list, tuple, set)):
                # 形式 1：{"Type": [2,5]}
                for cl in v:
                    mapping_cluster_to_type[str(cl)] = str(k)
            else:
                # 形式 2：{"2": "Type"}
                mapping_cluster_to_type[str(k)] = str(v)
        adata.obs[newcell_col] = leiden_as_str[subset_mask].map(mapping_cluster_to_type).astype("category")
    else:
        adata.obs[newcell_col] = fixed_label or "Adata1_Selected"
        adata.obs[newcell_col] = adata.obs[newcell_col].astype("category")

    return adata


def run_harmony_inplace(adata: ad.AnnData, vars_use: list[str], n_pcs: int) -> None:
    """在已有 PCA 基础上运行 Harmony，并写入 obsm['X_pca_harmony']。"""
    if not _HAS_HARMONY:
        raise RuntimeError("harmonypy not available. Please install it or set RUN_HARMONY=False.")
    if "X_pca" not in adata.obsm:
        raise ValueError("X_pca missing; run PCA first.")
    Z = adata.obsm["X_pca"]
    ho = hm.run_harmony(Z, adata.obs, vars_use, random_state=RANDOM_SEED)
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T


def hvg_scale_pca(adata: ad.AnnData, n_top: int, n_pcs: int) -> None:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, subset=False)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack")


def neighbors_cluster_umap(
    adata: ad.AnnData,
    n_pcs: int,
    n_neighbors: int,
    resolution: float,
    use_harmony: bool,
) -> None:
    if use_harmony and "X_pca_harmony" in adata.obsm:
        sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=RANDOM_SEED)
    else:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=RANDOM_SEED)
    sc.tl.leiden(adata, resolution=resolution, random_state=RANDOM_SEED)
    sc.tl.umap(adata, random_state=RANDOM_SEED)


def plot_umap_colors(adata: ad.AnnData, colors: list[str], prefix: str, dpi: int = 150, size: int = 3):
    # 仅绘制存在的 obs 列
    cols_exist = [c for c in colors if c in adata.obs]
    if not cols_exist:
        print("[UMAP] No valid obs columns found to color.")
        return
    sc.pl.umap(
        adata,
        color=cols_exist,
        ncols=min(3, len(cols_exist)),
        size=size,
        wspace=0.35,
        show=False,
        save=f"_{prefix}_umap.png",
        
    )


def plot_feature_umap(adata: ad.AnnData, genes: list[str], prefix: str, dpi: int = 150, size: int = 3):
    # 仅绘制存在于 var_names 的基因
    g_exist = [g for g in genes if g in adata.var_names]
    g_missing = [g for g in genes if g not in adata.var_names]
    if g_missing:
        print(f"[FeaturePlot] Missing genes (skipped): {g_missing}")
    if not g_exist:
        print("[FeaturePlot] No valid genes to plot.")
        return
    sc.pl.umap(
        adata,
        color=g_exist,
        ncols=4,
        size=size,
        wspace=0.35,
        vmax=None,
        show=False,
        save=f"_{prefix}_feature_umap.png",
        
    )


def cell_barplot(
    adata: ad.AnnData,
    groupby: str,
    hue: str = NEWCELL_COL,
    prefix: str = PREFIX,
    dpi: int = 150,
    percentage: bool = True,
):
    if hue not in adata.obs or groupby not in adata.obs:
        print(f"[Barplot] '{hue}' or '{groupby}' not in adata.obs; skip.")
        return
    df = adata.obs[[groupby, hue]].copy()
    ct = df.groupby([groupby, hue]).size().unstack(fill_value=0).sort_index()
    if percentage:
        ct = ct.div(ct.sum(axis=1), axis=0) * 100
        ylabel = "Percentage (%)"
    else:
        ylabel = "Cell count"

    fig, ax = plt.subplots(figsize=(10, 6), dpi=dpi)
    bottom = np.zeros(ct.shape[0])
    for col in ct.columns:
        ax.bar(ct.index.astype(str), ct[col].values, bottom=bottom, label=str(col))
        bottom += ct[col].values
    ax.set_title(f"{hue} composition by {groupby}")
    ax.set_xlabel(groupby)
    ax.set_ylabel(ylabel)
    ax.legend(title=hue, bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    out = os.path.join(OUTPUT_DIR, f"{prefix}_bar_{groupby}.png")
    fig.savefig(out)
    plt.close(fig)
    print(f"[Barplot] Saved: {out}")

def do_deg(adata, sample_prefix="sample",groupkey='leiden'):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby=groupkey, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=groupkey,
        n_genes=5,
        save=f'{sample_prefix}-holo_dotplot_leiden.png',
        min_logfoldchange=0.25
    )
    
    df1 = sc.get.rank_genes_groups_df(adata,group=None,log2fc_min=0.1,pval_cutoff=0.01)
    df1.to_csv(f'{sample_prefix}-holo_ranked_genes_{groupkey}.csv', index=False)
    print("DE analysis completed and results saved.")

# =========================
# Main
# =========================
if __name__ == "__main__":
    np.random.seed(RANDOM_SEED)
    _ensure_outdir(OUTPUT_DIR)

    print(f"[IO] Reading: {ADATA1_PATH}")
    adata1 = sc.read_h5ad(ADATA1_PATH)
    print(f"[IO] Reading: {ADATA2_PATH}")
    adata2 = sc.read_h5ad(ADATA2_PATH)

    # --- adata1：子集化并设置 newcelltype ---
    adata1 = subset_and_set_newcelltype(
        adata1,
        leiden_col=LEIDEN_COL,
        clusters_to_keep=CLUSTERS_TO_KEEP,
        newcelltype_mapping=NEWCELLTYPE_MAPPING,
        fixed_label=FIXED_NEWCELLTYPE_LABEL,
        newcell_col=NEWCELL_COL,
    )

    # --- adata2：检查 newcelltype 是否存在 ---
    if NEWCELL_COL not in adata2.obs:
        raise ValueError(f"'{NEWCELL_COL}' not found in adata2.obs. Please ensure adata2 already has this column.")

    # --- 在合并前确保两者都有 normalized 层 ---
    adata1 = ensure_normalized_layer(adata1, "normalized")
    adata2 = ensure_normalized_layer(adata2, "normalized")
    # --- 在合并前删除两份对象中的指定 obs 列 ---
    drop_obs_columns(adata1, OBS_COLS_TO_DROP)
    drop_obs_columns(adata2, OBS_COLS_TO_DROP)
    # --- 合并 ---
    print("[Merge] Concatenating adata1_subset and adata2 ...")
    merged = ad.concat([adata1, adata2], join="outer", label="batch", keys=["adata1", "adata2"], index_unique="-")

    # --- 用 normalized 层替换 X，并继续标准流程 ---
    merged = set_X_from_layer(merged, "normalized")

    # HVG & PCA
    hvg_scale_pca(merged, n_top=N_TOP_GENES, n_pcs=N_PCS)

    # Harmony（可选）
    use_harmony = RUN_HARMONY and _HAS_HARMONY and len(HARMONY_VARS) > 0 and all(v in merged.obs for v in HARMONY_VARS)
    if RUN_HARMONY and not _HAS_HARMONY:
        print("[Harmony] harmonypy not installed; skipping Harmony.")
    if RUN_HARMONY and _HAS_HARMONY and (not all(v in merged.obs for v in HARMONY_VARS)):
        print(f"[Harmony] Vars missing in obs: {[v for v in HARMONY_VARS if v not in merged.obs]} ; skipping Harmony.")
    if use_harmony:
        print(f"[Harmony] Running on vars: {HARMONY_VARS}")
        run_harmony_inplace(merged, vars_use=HARMONY_VARS, n_pcs=N_PCS)

    # Neighbors / Leiden / UMAP
    neighbors_cluster_umap(
        merged,
        n_pcs=N_PCS,
        n_neighbors=N_NEIGHBORS,
        resolution=LEIDEN_RES,
        use_harmony=use_harmony,
    )

    # --- UMAP（按指定列） ---
    plot_umap_colors(merged, [NEWCELL_COL, STAGE_COL, SAMPLE_COL], prefix=f"{PREFIX}_meta", dpi=FIG_DPI, size=POINT_SZ)

    # --- Feature UMAP ---
    plot_feature_umap(merged, FEATURE_GENES, prefix=f"{PREFIX}_feature", dpi=FIG_DPI, size=POINT_SZ)

    # --- Cell barplots（按 stage / sample） ---
    cell_barplot(merged, groupby=STAGE_COL, hue=NEWCELL_COL, prefix=PREFIX, dpi=FIG_DPI, percentage=True)
    cell_barplot(merged, groupby=SAMPLE_COL, hue=NEWCELL_COL, prefix=PREFIX, dpi=FIG_DPI, percentage=True)
    do_deg(merged, sample_prefix=dataset,groupkey='newcelltype')
    # --- 保存合并后的对象 ---
    out_h5ad = os.path.join(OUTPUT_DIR, MERGED_NAME)
    merged.write(out_h5ad)
    print(f"[Save] Merged AnnData saved to: {out_h5ad}")

    print("[Done] All tasks completed.")
