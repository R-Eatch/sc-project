#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
samap_weighted_trend_pipeline.py  (WEIGHTED ONLY)
-------------------------------------------------
功能：
1) 读取 AnnData (.h5ad)，用 adata.obs['celltype'] 计算各细胞类型比例（权重），保存 CSV。
2) 读取三张 SAMap 细胞类型×细胞类型矩阵（msMG→sgMG / sgAG / sgSG）。
3) 将三腺体目标列合并，按“同名 celltype”计算每个源类型的：
   - 点双列相关 r_pb（scores vs gland label）。
4) 仅用“加权”做统计检验：
   - 加权 Page’s 趋势检验（SG < AG < MG）
   - 计划性对比（MG>AG, AG>SG）的加权一侧置换检验（sign-flip）
5) 可视化：
   - r_pb 的箱线图（SG/AG/MG），在图上画显著性括号、星号、p 值，并在图注写明方法与关键数值
   - r_pb 趋势“意大利面”图（叠加加权均值曲线）
6) 保存所有结果表与图片。

依赖：
    pip install anndata pandas numpy matplotlib
（scipy 可选，不依赖）
"""

# ============ 全局变量（按需改这里） ============
ADATA_PATH = "../M-MG_counts_pr.h5ad"   # path to your AnnData (.h5ad)
OBS_CELLTYPE_COL = "celltype"           # the obs column to use for weights

SAMAP_SG_MG_PATH = "../celltype/msMG-sgMG-celltype-MappingTable.csv"
SAMAP_SG_AG_PATH = "../celltype/msMG-sgAG-celltype-MappingTable.csv"
SAMAP_SG_SG_PATH = "../celltype/msMG-sgSG-celltype-MappingTable.csv"

OUTPUT_PREFIX = "SG_corr_WEIGHTED"      # 输出文件前缀
N_PERM_SIGNFLIP = 50000                 # 置换次数（配对差）
N_PERM_PAGE = 20000                     # 置换次数（Page 趋势）
RANDOM_SEED = 7                         # 随机种子
# ============================================

from typing import Dict, Tuple, Optional
import unicodedata, re, os, json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- 基础工具 ----------
def norm_name(x: str) -> str:
    s = str(x)
    s = unicodedata.normalize("NFKC", s).strip()
    s = re.sub(r"\s+", " ", s)
    return s

def strip_prefix(x: str, pref: str) -> str:
    s = norm_name(x)
    if s.lower().startswith(pref.lower()):
        return s[len(pref):].lstrip("_ ").strip()
    return s

def p_to_stars(p: float) -> str:
    if p is None or not np.isfinite(p): return "n.s."
    if p < 1e-4: return "****"
    if p < 1e-3: return "***"
    if p < 1e-2: return "**"
    if p < 0.05: return "*"
    return "n.s."

# ---------- 读数据 ----------
def read_h5ad_weights(adata_path: str, obs_col: str) -> pd.DataFrame:
    try:
        import anndata as ad
    except Exception as e:
        raise RuntimeError("需要 anndata：pip install anndata") from e
    if not os.path.exists(adata_path):
        raise FileNotFoundError(f"h5ad 不存在: {adata_path}")
    adata = ad.read_h5ad(adata_path)
    if obs_col not in adata.obs.columns:
        raise KeyError(f"'{obs_col}' 不在 adata.obs 中。可用列：{list(adata.obs.columns)}")
    ct = adata.obs[obs_col].astype(str).map(norm_name)
    counts = ct.value_counts(dropna=False).rename_axis("celltype").to_frame("count").reset_index()
    counts["proportion"] = counts["count"] / counts["count"].sum()
    counts = counts.sort_values("proportion", ascending=False).reset_index(drop=True)
    return counts

def read_samap_matrix(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"SAMap 矩阵不存在: {path}")
    df = pd.read_csv(path, sep=None, engine="python", index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df.index = [norm_name(i) for i in df.index]
    df.columns = [norm_name(c) for c in df.columns]
    return df

# ---------- 构建合并目标 + 计算指标 ----------
def build_combined_target_and_rpb(MG: pd.DataFrame, AG: pd.DataFrame, SG: pd.DataFrame
                                  ) -> Tuple[pd.DataFrame, pd.DataFrame]:
    common_rows = sorted(set(MG.index) & set(AG.index) & set(SG.index))
    if not common_rows:
        raise RuntimeError("三张 SAMap 矩阵没有共同的源行。")
    MGc = MG.loc[common_rows].copy(); MGc.columns = [f"MG::{c}" for c in MGc.columns]
    AGc = AG.loc[common_rows].copy(); AGc.columns = [f"AG::{c}" for c in AGc.columns]
    SGc = SG.loc[common_rows].copy(); SGc.columns = [f"SG::{c}" for c in SGc.columns]
    combined = pd.concat([SGc, AGc, MGc], axis=1)

    gland_for_col = pd.Series(index=combined.columns, dtype=object)
    for c in combined.columns:
        gland_for_col[c] = c.split("::", 1)[0]  # SG / AG / MG

    def point_biserial(scores: np.ndarray, labels01: np.ndarray) -> float:
        mask = np.isfinite(scores) & np.isfinite(labels01)
        x = scores[mask]; y = labels01[mask].astype(float)
        if x.size < 3: return np.nan
        if y.min() == y.max(): return np.nan
        x_mean = x.mean(); y_mean = y.mean()
        x_std = x.std(ddof=1); y_std = y.std(ddof=1)
        if x_std == 0 or y_std == 0: return np.nan
        cov = ((x - x_mean) * (y - y_mean)).sum() / (x.size - 1)
        return float(cov / (x_std * y_std))

    rows = []
    for r in common_rows:
        s = combined.loc[r].to_numpy(dtype=float)
        r_pb = {}
        for g in ("SG","AG","MG"):
            y = (gland_for_col.values == g).astype(int)
            r_pb[g] = point_biserial(s, y)
        rows.append({"source_celltype": r,
                     "r_pb_SG": r_pb["SG"], "r_pb_AG": r_pb["AG"], "r_pb_MG": r_pb["MG"]})
    rpb = pd.DataFrame(rows).set_index("source_celltype")
    return rpb, combined

def attach_weights_by_basename(rpb: pd.DataFrame, weights: pd.DataFrame) -> pd.DataFrame:
    wmap = dict(zip(weights["celltype"].map(norm_name), weights["proportion"]))
    base = rpb.index.to_series().map(lambda s: strip_prefix(s, "ms_"))
    w = base.map(lambda b: wmap.get(norm_name(b), np.nan))
    out = rpb.copy()
    out.insert(0, "base_name", base.values)
    out.insert(1, "weight", w.values)
    return out

# ---------- 统计检验（加权） ----------
def signflip_perm_weighted_mean(d: np.ndarray, w: np.ndarray, iters: int, seed: int) -> Tuple[float, float]:
    mask = np.isfinite(d) & np.isfinite(w)
    d = d[mask]; w = w[mask]
    if d.size == 0: return np.nan, np.nan
    w = w / (w.sum() if w.sum()!=0 else 1.0)
    obs = float(np.sum(w * d))
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1,1], size=(iters, d.size))
    sim = (w * d)[None, :] * signs
    p = float((sim.sum(axis=1) >= obs).mean())
    return obs, p

def pages_L_weighted(mat: np.ndarray, w: np.ndarray, iters: int, seed: int) -> Tuple[float, float]:
    """加权 Page’s L，mat 形状 (n,3) 按 [SG, AG, MG]。"""
    mask_rows = np.isfinite(mat).all(axis=1) & np.isfinite(w)
    X = mat[mask_rows]; w = w[mask_rows]
    if X.shape[0] == 0: return np.nan, np.nan
    ranks = np.apply_along_axis(lambda v: pd.Series(v).rank(method="average").to_numpy(), 1, X)
    pos_w = np.arange(1,4, dtype=float)  # 1,2,3
    L_obs = float(np.sum(w[:,None] * ranks * pos_w))
    perms = np.array([[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]], dtype=int)
    row_perm_scores = np.empty((X.shape[0], 6), dtype=float)
    for i in range(X.shape[0]):
        v = X[i]
        for p_idx, perm in enumerate(perms):
            r = pd.Series(v[perm]).rank(method="average").to_numpy()
            row_perm_scores[i, p_idx] = float((r * pos_w).sum())
    rng = np.random.default_rng(seed)
    choice_idx = rng.integers(0, 6, size=(iters, X.shape[0]))  # (iters, n)
    sim_vals = (row_perm_scores[np.arange(X.shape[0])[:,None], choice_idx.T] * w[:,None]).sum(axis=0)
    p = float((sim_vals >= L_obs).mean())
    return L_obs, p

# ---------- 画图（含显著性标注） ----------
def annotate_sig(ax, x1, x2, y, text, h=0.02):
    """在箱线图上方画显著性括号与文本。x1/x2 是箱体索引（0基），y 为起始 y 值。"""
    ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.2, c='k')
    ax.text((x1+x2)/2, y+h*1.2, text, ha='center', va='bottom')

def save_boxplot_with_significance(rpb_w: pd.DataFrame, summary: Dict, out_png: str):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    data = [rpb_w["r_pb_SG"].to_numpy(), rpb_w["r_pb_AG"].to_numpy(), rpb_w["r_pb_MG"].to_numpy()]
    bp = ax.boxplot(data, labels=["SG","AG","MG"], showfliers=False)
    ax.set_ylabel("Point-biserial correlation (r_pb)")
    ax.set_title("Weighted analysis (per cell type)")

    # 取三个组的最大 y 用于放括号
    ymax = np.nanmax([np.nanmax(d) for d in data])
    ymin = np.nanmin([np.nanmin(d) for d in data])
    span = (ymax - ymin) if np.isfinite(ymax - ymin) and (ymax - ymin) > 0 else 0.05
    base_y = ymax + span*0.1
    step = span*0.12

    # 标注两两对比（加权）的 p 与星号
    p_ag_sg = summary["paired_diffs_weighted"]["AG-SG"]["p_perm_one_sided"]
    p_mg_ag = summary["paired_diffs_weighted"]["MG-AG"]["p_perm_one_sided"]
    annotate_sig(ax, 0+1-1, 1+1-1, base_y, f"AG > SG: {p_to_stars(p_ag_sg)}  (p={p_ag_sg:.3g})")
    annotate_sig(ax, 1+1-1, 2+1-1, base_y+step, f"MG > AG: {p_to_stars(p_mg_ag)}  (p={p_mg_ag:.3g})")

    # 在图底部写图注简述（方法 + Page’s p + 加权均值）
    means = summary["r_pb_group_means_weighted"]
    p_page = summary["page_L_weighted"]["p_one_sided"]
    caption = (f"Weighted Page’s trend test for SG<AG<MG: p={p_page:.3g}.  "
               f"Weighted means r_pb: SG={means['SG']:.4f}, AG={means['AG']:.4f}, MG={means['MG']:.4f}.  "
               f"Paired sign-flip tests are one-sided on weighted means.")
    fig.text(0.02, 0.02, caption, ha='left', va='bottom', fontsize=9)
    plt.tight_layout(rect=[0,0.06,1,1])
    plt.savefig(out_png, dpi=200)
    plt.close()

def save_spaghetti_rpb(rpb_w: pd.DataFrame, out_png: str, weights: np.ndarray):
    xs = [0,1,2]
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    for _, row in rpb_w[["r_pb_SG","r_pb_AG","r_pb_MG"]].iterrows():
        y = row.to_numpy(dtype=float)
        ax.plot(xs, y, alpha=0.6)
    # 加权均值线
    vals = rpb_w[["r_pb_SG","r_pb_AG","r_pb_MG"]].to_numpy(float)
    m = np.isfinite(vals).all(axis=1) & np.isfinite(weights)
    w = weights[m]
    w = w / (w.sum() if w.sum()!=0 else 1.0)
    sub = rpb_w.loc[m]
    mean_line = np.array([np.sum(w * sub["r_pb_SG"].to_numpy()),
                          np.sum(w * sub["r_pb_AG"].to_numpy()),
                          np.sum(w * sub["r_pb_MG"].to_numpy())])
    ax.plot(xs, mean_line, lw=3)
    ax.set_xticks(xs); ax.set_xticklabels(["SG","AG","MG"])
    ax.set_ylabel("Point-biserial correlation (r_pb)")
    ax.set_title("Per-source r_pb trend: SG → AG → MG (weighted mean line)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

# ---------- 主流程 ----------
def run_pipeline():
    # 1) 权重
    weights = read_h5ad_weights(ADATA_PATH, OBS_CELLTYPE_COL)
    w_csv = f"{OUTPUT_PREFIX}.weights_from_{OBS_CELLTYPE_COL}.csv"
    weights.to_csv(w_csv, index=False)

    # 2) SAMap 矩阵
    MG = read_samap_matrix(SAMAP_SG_MG_PATH)
    AG = read_samap_matrix(SAMAP_SG_AG_PATH)
    SG = read_samap_matrix(SAMAP_SG_SG_PATH)

    # 3) 合并目标 + r_pb
    rpb, combined = build_combined_target_and_rpb(MG, AG, SG)
    rpb_csv = f"{OUTPUT_PREFIX}.rpb_per_row.csv"; rpb.to_csv(rpb_csv)

    # 4) 同名映射附加权重
    rpb_w = attach_weights_by_basename(rpb, weights)
    rpb_w_csv = f"{OUTPUT_PREFIX}.rpb_with_weights.csv"; rpb_w.to_csv(rpb_w_csv)

    # 5) 统计（仅加权）
    mat = rpb_w[["r_pb_SG","r_pb_AG","r_pb_MG"]].to_numpy(float)
    w_vec = rpb_w["weight"].to_numpy(float)

    # Page’s 趋势（加权）
    L_w, pL_w = pages_L_weighted(mat, w_vec, iters=N_PERM_PAGE, seed=RANDOM_SEED+3)

    # 计划性对比（加权）
    d_MG_AG = rpb_w["r_pb_MG"].to_numpy(float) - rpb_w["r_pb_AG"].to_numpy(float)
    d_AG_SG = rpb_w["r_pb_AG"].to_numpy(float) - rpb_w["r_pb_SG"].to_numpy(float)
    obs_w_1, p_w_1 = signflip_perm_weighted_mean(d_MG_AG, w_vec, iters=N_PERM_SIGNFLIP, seed=RANDOM_SEED+11)
    obs_w_2, p_w_2 = signflip_perm_weighted_mean(d_AG_SG, w_vec, iters=N_PERM_SIGNFLIP, seed=RANDOM_SEED+13)

    # 组加权均值
    def wmean(vals, w):
        m = np.isfinite(vals) & np.isfinite(w)
        if m.sum()==0: return float("nan")
        ww = w[m]; vv = vals[m]; s = ww.sum()
        return float(np.sum(ww*vv)/s) if s>0 else float("nan")

    means_w = {
        "SG": wmean(rpb_w["r_pb_SG"].to_numpy(float), w_vec),
        "AG": wmean(rpb_w["r_pb_AG"].to_numpy(float), w_vec),
        "MG": wmean(rpb_w["r_pb_MG"].to_numpy(float), w_vec),
    }

    # 汇总 JSON
    summary = {
        "n_rows_total": int(rpb_w.shape[0]),
        "weights_csv": w_csv,
        "rpb_per_row_csv": rpb_csv,
        "rpb_with_weights_csv": rpb_w_csv,
        "r_pb_group_means_weighted": means_w,
        "page_L_weighted": {"L": L_w, "p_one_sided": pL_w},
        "paired_diffs_weighted": {
            "MG-AG": {"mean": float(obs_w_1), "p_perm_one_sided": float(p_w_1)},
            "AG-SG": {"mean": float(obs_w_2), "p_perm_one_sided": float(p_w_2)},
        },
        "notes": "All tests are one-sided, weighted by source cell-type proportions from the h5ad."
    }
    with open(f"{OUTPUT_PREFIX}.summary.json","w",encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    # 6) 画图（含显著性）
    save_boxplot_with_significance(rpb_w, summary, f"{OUTPUT_PREFIX}.rpb_boxplot_weighted.pdf")
    save_spaghetti_rpb(rpb_w, f"{OUTPUT_PREFIX}.rpb_spaghetti_weighted.pdf", w_vec)

    print("[DONE] Saved outputs with prefix:", OUTPUT_PREFIX)
    print("Summary JSON:", f"{OUTPUT_PREFIX}.summary.json")
    print("Weights CSV:", w_csv)
    print("r_pb per-row CSV:", rpb_csv)
    print("r_pb with weights CSV:", rpb_w_csv)
    print("Figures:",
          f"{OUTPUT_PREFIX}.rpb_boxplot_weighted.png",
          f"{OUTPUT_PREFIX}.rpb_spaghetti_weighted.png")

if __name__ == "__main__":
    run_pipeline()

