#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Compute cross-gland mean expression for a *celltype-aware* conserved-gene list,
then build a unique Top-N `selection` dict (one gene belongs to one celltype only).

Inputs
------
1. unique_common_genes_long.csv  ← long-format, columns: celltype,gene,[...]
2. single_adata_normalized.h5ad  ← ONE file, with .X as normalized data
   └── AnnData.obs[GROUP_KEY]     ← cell-type annotation (e.g., 'newcelltype')
   └── AnnData.obs[GLAND_KEY]     ← gland annotation (e.g., 'gland')

Outputs
-------
1. cross_gland_gene_expression.csv ← long table: celltype,gene,cross_gland_mean,EG_mean,MG_mean,SG_mean
2. Prints + returns `selection`      ← {celltype: [Top-N genes]}
"""

from pathlib import Path
from collections import OrderedDict
import numpy as np
import pandas as pd
import scanpy as sc

# --------------------------- CONFIG ---------------------------
### (已修改)
DATA_ROOT   = Path(r"./")
# ↓↓↓ 你的单个 .h5ad 文件名
H5AD_FILE   = "../../h5ad2seurat/processed/v9_all_log1p.h5ad" 
# ↓↓↓ 你的基因列表 (来自上一步)
GENE_CSV    = Path("unique_common_genes_long.csv") 

# ↓↓↓ 预处理后你关心的腺体 (必须与 rename_and_merge_glands 的输出匹配)
GLANDS_TO_ANALYZE = ["EG", "MG", "SG"] 

# ↓↓↓ AnnData.obs 中用于细胞类型和腺体的列名
GROUP_KEY   = "anno"   # 细胞类型 (e.g., MaSC, Luminal)
GLAND_KEY   = "gland"         # 腺体 (e.g., MG, SG, AG, CG)

TOP_N       = 5
CSV_OUT     = "cross_gland_gene_expression.csv" # (已重命名)
# ---------------------------------------------------------------


# ---------- 1. (新增) 你的预处理函数 ----------
def rename_and_merge_glands(adata):
    """Apply:
     1) SG -> EG
     2) AG -> SG, CG -> SG  (after step 1)
    """
    if GLAND_KEY not in adata.obs.columns:
        raise ValueError(f"缺少 adata.obs['{GLAND_KEY}'] 列。")

    print(f"[INFO] 正在重命名/合并 '{GLAND_KEY}' 列...")
    print(f"  原始计数:\n{adata.obs[GLAND_KEY].value_counts().to_string()}")

    # Step 1: SG -> EG
    adata.obs[GLAND_KEY] = adata.obs[GLAND_KEY].replace({'SG': 'EG'})
    # Step 2: merge AG/CG into SG
    adata.obs[GLAND_KEY] = adata.obs[GLAND_KEY].replace({'AG': 'SG', 'CG': 'SG'})

    print("[OK] 已重命名/合并 gland：SG→EG；AG+CG→SG。")
    print(f"[INFO] 重命名后 '{GLAND_KEY}' 计数：")
    print(adata.obs[GLAND_KEY].value_counts().to_string())
    return adata


# ---------- 2. (不变) 读取 celltype↔gene 映射 ----------
def load_mapping(csv_path: Path) -> dict[str, list[str]]:
    """
    读取CSV，构建一个 {celltype -> [genes]} 的字典。
    关键：通过'seen'集合确保每个基因只被分配给它在CSV中*首次*出现的celltype。
    """
    df = pd.read_csv(csv_path)

    df.columns = [str(c).lower() for c in df.columns]

    # 自动检测 'celltype' 和 'gene' 列，否则使用前两列
    col_map = {}
    if "celltype" in df.columns:
        col_map["celltype"] = "celltype"
    else:
        print(f"[WARN] 未找到 'celltype' 列，使用第1列 '{df.columns[0]}'")
        col_map["celltype"] = df.columns[0]
        
    if "gene" in df.columns:
        col_map["gene"] = "gene"
    else:
        print(f"[WARN] 未找到 'gene' 列，使用第2列 '{df.columns[1]}'")
        col_map["gene"] = df.columns[1]

    if df.shape[1] < 2:
        raise ValueError("CSV needs at least two columns (celltype, gene).")
    df = df.rename(columns={v: k for k, v in col_map.items()})


    mapping: dict[str, list[str]] = OrderedDict()
    seen = set()
    for ct, g in zip(df["celltype"], df["gene"]):
        if pd.isna(ct) or pd.isna(g) or g in seen:
            continue
        mapping.setdefault(ct, []).append(g)
        seen.add(g)
    return mapping


# ---------- 3. (已修改) 计算指定细胞类型/腺体中一组基因的均值 ----------
def mean_expr_gland(adata: sc.AnnData, genes: list[str], celltype: str, gland: str) -> pd.Series:
    """
    (替换了 mean_expr_species)
    返回一个物种中，一个 celltype 在一个 gland 中的基因平均表达。
    """
    
    # 检查细胞类型和腺体是否存在于 .obs 中
    if celltype not in adata.obs[GROUP_KEY].values:
        print(f"  [WARN] {celltype} 不在 {GROUP_KEY} 中。")
        return pd.Series(np.nan, index=genes)
    if gland not in adata.obs[GLAND_KEY].values:
        print(f"  [WARN] {gland} 不在 {GLAND_KEY} 中。")
        return pd.Series(np.nan, index=genes)

    # 取目标细胞 (必须同时匹配 celltype 和 gland)
    idx_cells = (adata.obs[GROUP_KEY] == celltype) & (adata.obs[GLAND_KEY] == gland)
    
    if idx_cells.sum() == 0:
        # print(f"  [Debug] {celltype} + {gland} 中没有细胞。")
        return pd.Series(np.nan, index=genes)

    # 限制到存在于该 AnnData 的基因，保持原顺序
    valid_genes = [g for g in genes if g in adata.var_names]
    if not valid_genes:
        print(f"  [WARN] {celltype} 的基因列表在 anndata.var_names 中均未找到。")
        return pd.Series(np.nan, index=genes)

    # (已修改) 从 .X 提取表达矩阵并求均值
    X = adata[idx_cells, valid_genes].X  # 使用 .X
    
    means = np.asarray(X.mean(axis=0)).ravel()
    return pd.Series(means, index=valid_genes).reindex(genes) # 保证长度一致


# ---------- 4. (已修改) 构建长格式表达表 ----------
def build_expression_table(adata: sc.AnnData, mapping: dict) -> pd.DataFrame:
    """
    (替换了 adatas: dict)
    遍历所有 celltype，计算其在每个 gland 中的平均表达，然后计算跨腺体均值。
    """
    recs = []
    print("\n[INFO] 正在计算所有 celltype 和 gland 组合的平均表达...")
    
    for ct, genes in mapping.items():
        # print(f"  正在处理 {ct}...") # (如果需要详细日志，取消注释)
        
        # (已修改) 不再遍历 adatas 字典，而是遍历 GLANDS_TO_ANALYZE 列表
        gland_means = {
            gland: mean_expr_gland(adata, genes, ct, gland)
            for gland in GLANDS_TO_ANALYZE
        }
        
        stacked = pd.DataFrame(gland_means)
        
        # (已修改) 计算 cross_gland_mean
        stacked["cross_gland_mean"] = stacked.mean(axis=1, skipna=True)
        
        # (已修改) 动态写入记录
        recs.extend(
            {
                "celltype": ct,
                "gene": g,
                "cross_gland_mean": row["cross_gland_mean"],
                # 动态添加每个gland的均值
                **{f"{gland}_mean": row.get(gland, np.nan) for gland in GLANDS_TO_ANALYZE}
            }
            for g, row in stacked.iterrows()
        )
    print("[OK] 表达表计算完成。")
    return pd.DataFrame.from_records(recs)


# ---------- 5. (已修改) 生成 Top-N selection ----------
def build_selection(expr_df: pd.DataFrame, mapping: dict) -> dict:
    """
    基于 'cross_gland_mean' 进行排序和筛选。
    """
    selection = {}
    used = set()  # 已入选基因
    print("\n[INFO] 正在构建唯一的 Top-N 基因选择...")
    
    # 保持 mapping 中 celltype 的顺序 (优先级)
    for ct in mapping:
        top_genes = (
            expr_df.query("celltype == @ct and gene not in @used")
                   # (已修改) 按 cross_gland_mean 排序
                   .sort_values("cross_gland_mean", ascending=False)
                   .head(TOP_N)["gene"]
                   .tolist()
        )
        selection[ct] = top_genes
        used.update(top_genes)
        
    print("[OK] Top-N 字典构建完成。")
    return selection


# ---------- 6. (已修改) 主流程 ----------
def main():
    # 1) 读取映射
    mapping_path = DATA_ROOT / GENE_CSV
    mapping = load_mapping(mapping_path)
    print(f"[OK] 基因映射加载完成: {sum(map(len, mapping.values()))} 个唯一基因 "
          f"分布在 {len(mapping)} 个 celltype 中 (来自 {mapping_path.name})")

    # 2) 载入单个 AnnData
    adata_path = DATA_ROOT / H5AD_FILE
    print(f"[INFO] 正在加载 AnnData: {adata_path} ...")
    ad = sc.read_h5ad(adata_path)
    print(f"[OK] AnnData 加载: {ad.n_obs:,} 个细胞, {ad.n_vars:,} 个基因。")
    
    # 验证 .obs 列
    if GROUP_KEY not in ad.obs.columns:
        raise KeyError(f"{adata_path} - 缺少细胞类型列 .obs['{GROUP_KEY}']")
    if GLAND_KEY not in ad.obs.columns:
        raise KeyError(f"{adata_path} - 缺少腺体列 .obs['{GLAND_KEY}']")
        
    print(f"[INFO] 将使用 adata.X 作为归一化表达数据。")

    # 3) (新增) 运行预处理
    ad = rename_and_merge_glands(ad)

    # 4) 表达量统计
    expr_df = build_expression_table(ad, mapping)
    
    output_path = DATA_ROOT / CSV_OUT
    expr_df.to_csv(output_path, index=False)
    print(f"\n[OK] 跨腺体表达CSV已保存 → {output_path} ({len(expr_df):,} 行)")

    # 5) Top-N selection
    selection = build_selection(expr_df, mapping)
    print("\n" + "="*30)
    print(f"    Top-{TOP_N} 'selection' 字典    ")
    print("="*30)
    for ct, genes in selection.items():
        print(f"{ct:<15} : {genes}")

    return selection


if __name__ == "__main__":
    # 确保脚本在DATA_ROOT目录运行，或者GENE_CSV和H5AD_FILE在正确路径
    try:
        selection = main()
    except FileNotFoundError as e:
        print(f"\n[ERROR] 文件未找到: {e}")
        print("请检查 CONFIG 部分中的 DATA_ROOT, H5AD_FILE, 和 GENE_CSV 路径是否正确。")
    except KeyError as e:
        print(f"\n[ERROR] 键错误: {e}")
        print(f"请检查 CONFIG 中的 GROUP_KEY ('{GROUP_KEY}') 或 GLAND_KEY ('{GLAND_KEY}') "
              "是否与 .h5ad 文件中的 .obs 列名匹配。")
