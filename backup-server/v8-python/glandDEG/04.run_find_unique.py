#!/usr/bin/env python
# coding: utf-8
"""

"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# --- (已修改) ---
# 你的单个 AnnData 文件
H5AD_FILE = Path("../../h5ad2seurat/processed/v9_all_log1p.h5ad")

# 要分析的腺体 (在预处理后)
GLANDS = ["EG", "MG", "SG"]

# 对应 DEG csv 路径 (!! 假设你已按腺体准备了这些文件 !!)
DEG_PATH = {
    gland: Path(f"deg_out/MyDataset_DEG_{gland}.csv") 
    for gland in GLANDS
}
# --- (结束修改) ---

GROUP_KEY  = "anno" # 细胞类型列
GLAND_KEY  = "gland"       # 腺体列 (拆分会用到)
ZERO_EXPR_THRESHOLD = 0.01
TOP_N      = 5
HEAT_CMAP  = "flare_r"
DPI        = 300

# Gland 调色板
GLAND_CMAP = {
    "EG": sns.light_palette("#7BC6E6", n_colors=256, as_cmap=True), # 蓝色
    "MG": sns.light_palette("#E5B764", n_colors=256, as_cmap=True), # 黄色
    "SG": sns.light_palette("#A3C29F", n_colors=256, as_cmap=True), # 绿色
}
# ------------------------------------------------------------------


# 腺体预处理函数
def rename_and_merge_glands(adata, gland_col):
    """
    应用: 1) SG -> EG; 2) AG -> SG, CG -> SG
    """
    if gland_col not in adata.obs.columns:
        raise ValueError(f"缺少 adata.obs['{gland_col}'] 列。")
    
    print(f"[INFO] 正在重命名/合并 '{gland_col}' 列...")
    if pd.api.types.is_categorical_dtype(adata.obs[gland_col]):
        adata.obs[gland_col] = adata.obs[gland_col].astype(str)

    adata.obs[gland_col] = adata.obs[gland_col].replace({'SG': 'EG'})
    adata.obs[gland_col] = adata.obs[gland_col].replace({'AG': 'SG', 'CG': 'SG'})

    adata.obs[gland_col] = adata.obs[gland_col].astype("category")
    print(f"[OK] 腺体合并完成。")
    print(f"[INFO] 重命名后 '{gland_col}' 计数：\n{adata.obs[gland_col].value_counts().to_string()}")
    return adata


# (已修改) 加载、预处理并按 Gland 拆分
def load_and_split_anndata():
    """
    读 1 份 AnnData，预处理，按 Gland 拆分到 dict{gland: AnnData}
    (已移除 MaSC 合并)
    (新增) 返回自动检测到的细胞类型列表
    """
    print(f"[INFO] 正在加载主 AnnData: {H5AD_FILE} ...")
    f = H5AD_FILE
    if not f.exists():
        raise FileNotFoundError(f)
    
    ad_main = sc.read_h5ad(f)
    print(f"[OK] 主 AnnData 加载: {ad_main.n_obs:,} 个细胞, {ad_main.n_vars:,} 个基因。")

    # 1. (已移除) MaSC 合并逻辑
    
    # (新增) 自动检测所有细胞类型，并排序
    if pd.api.types.is_categorical_dtype(ad_main.obs[GROUP_KEY]):
        all_celltypes = ad_main.obs[GROUP_KEY].cat.categories.tolist()
    else:
        all_celltypes = sorted(ad_main.obs[GROUP_KEY].unique().tolist())
    print(f"[INFO] 自动检测到 {len(all_celltypes)} 个细胞类型 (按此排序):")
    print(all_celltypes)

    # 2. 运行腺体预处理
    ad_main = rename_and_merge_glands(ad_main, GLAND_KEY)
    

    # 4. 按 Gland 拆分
    ads_dict = {}
    for gland in GLANDS:
        adata_subset = ad_main[ad_main.obs[GLAND_KEY] == gland].copy()
        if adata_subset.n_obs > 0:
            ads_dict[gland] = adata_subset
            print(f"  > 拆分 '{gland}': {adata_subset.n_obs:,} 个细胞。")
        else:
            print(f"  [WARN] '{gland}' 中没有细胞，将跳过此组。")
            
    # (修改) 返回字典和细胞类型列表
    return ads_dict, all_celltypes


def genes_zero_in_others(genes, other_adatas):
    """(逻辑不变) 过滤：在所有“其他腺体”里 ≤ N % 细胞表达的基因集合."""
    keep = set(genes)
    for ad in other_adatas:
        genes_to_check = [g for g in keep if g in ad.var_names]
        if not genes_to_check:
            continue
            
        mat_other = (ad[:, genes_to_check].X > 0)
        prop = mat_other.mean(axis=0).ravel()
        gene_prop_map = dict(zip(genes_to_check, prop))
        
        keep -= {g for g in keep if gene_prop_map.get(g, 0) > ZERO_EXPR_THRESHOLD}
        
    return keep


# (已修改) 
def filter_gland_specific(gland, ads_dict, celltypes_to_keep):
    """
    (替换 filter_species_specific)
    (修改) 使用传入的 celltypes_to_keep 列表进行筛选
    """
    if gland not in ads_dict:
        print(f"  [WARN] {gland} 不在已加载的 anndata 字典中，跳过。")
        return {}
        
    ad_target = ads_dict[gland]
    others    = [ads_dict[g] for g in ads_dict if g != gland]

    zero_in_others = genes_zero_in_others(ad_target.var_names, others)
    print(f"    > {len(zero_in_others)} 个基因在其他腺体中低表达。")

    deg_file = DEG_PATH[gland]
    if not deg_file.exists():
        print(f"  [ERROR] 找不到DEG文件: {deg_file}。跳过 {gland}。")
        return {}
        
    deg = pd.read_csv(deg_file)
    deg = deg[deg["names"].isin(zero_in_others)]
    deg = deg.sort_values("logfoldchanges", ascending=False)

    # (修改) 按自动检测到的细胞类型列表筛选
    out = {ct: df for ct, df in deg.groupby("group") if ct in celltypes_to_keep}
    return out


def save_deg_csv(ct_dict, gland):
    """(逻辑不变)"""
    if not ct_dict:
        return None
    df = pd.concat(
        [d.assign(celltype=ct) for ct, d in ct_dict.items()],
        ignore_index=True
    )
    if df.empty:
        return None
    out = f"{gland}_specific_DEG.csv"
    df.to_csv(out, index=False)
    return out


# (已修改)
def top_genes(ct_dict, celltype_order_list, n=TOP_N):
    """(修改) 使用传入的 celltype_order_list 保证顺序"""
    genes = []
    # (修改) 按自动检测的顺序遍历
    for ct in celltype_order_list: 
        if ct in ct_dict:
            genes += ct_dict[ct].nlargest(n, 'logfoldchanges')['names'].tolist()
    
    seen = set(); ordered = []
    for g in genes:
        if g not in seen:
            ordered.append(g); seen.add(g)
    return ordered


# (已修改)
def plot_heatmap(genes, adata, gland_tag, gene_list_tag, celltype_order_list):
    """
    (修改) 使用传入的 celltype_order_list 保证行顺序
    """
    
    genes_avail = [g for g in genes if g in adata.var_names]
    if not genes_avail:
        print(f"  [WARN] 基因列表 {gene_list_tag} 在 {gland_tag} 数据中均不存在。")
        return None

    # (修改) 按自动检测的顺序筛选
    present_ct = [ct for ct in celltype_order_list if ct in adata.obs[GROUP_KEY].unique()]
    if not present_ct:
        print(f"  [WARN] {gland_tag} 数据中不存在 {celltype_order_list} 中的任何细胞类型。")
        return None

    expr = pd.DataFrame(index=genes_avail, columns=present_ct, dtype=float)
    var_lookup = pd.Series(np.arange(adata.n_vars), index=adata.var_names)

    for ct in present_ct:
        cells = adata.obs_names[adata.obs[GROUP_KEY] == ct]
        mat   = adata[cells].X
        mean_vec = np.asarray(mat.mean(axis=0)).ravel()
        expr[ct] = mean_vec[var_lookup[genes_avail].values]

    # Z-score (每基因行) 并截断到 [0,2]
    expr_z = (
        expr.sub(expr.mean(axis=1), axis=0)
            .div(expr.std(axis=1).replace(0, 1), axis=0)
            .clip(lower=0, upper=2) 
            .T
    ).loc[present_ct]               # (修改) 保证行顺序
    
    row_colors = None
    if f"{GROUP_KEY}_colors" in adata.uns:
        try:
            # 确保 categories 和 colors 列表能匹配
            if pd.api.types.is_categorical_dtype(adata.obs[GROUP_KEY]):
                cat = adata.obs[GROUP_KEY].cat.categories
            else:
                cat = sorted(adata.obs[GROUP_KEY].unique())
                
            palette  = adata.uns[f"{GROUP_KEY}_colors"]
            
            # 如果调色板长度与类别数不匹配，进行适配
            if len(palette) == len(cat):
                lut = dict(zip(cat, palette))
                row_colors = [lut.get(ct, '#808080') for ct in present_ct] # 使用 .get
            else:
                 print(f"  [WARN] 颜色数量 ({len(palette)}) 与类别数 ({len(cat)}) 不匹配。")

        except Exception as e:
            print(f"  [WARN] 提取行颜色失败: {e}")

    # --- 选 Gland 专属色板 (基于基因来源) ---
    cmap_choice = GLAND_CMAP.get(gene_list_tag, HEAT_CMAP)

    sns.set(style="white", font_scale=0.75)
    g = sns.clustermap(
        expr_z,
        row_cluster=False, col_cluster=False,
        cmap=cmap_choice,
        vmin=0, vmax=2,           
        row_colors=row_colors,
        linewidths=0.1,
        figsize=(max(6, len(genes_avail)*0.4), max(3, len(present_ct)*0.45))
    )
    g.ax_heatmap.set_xlabel("Gene")
    g.ax_heatmap.set_ylabel("Cell type")
    g.ax_heatmap.set_title(f"{gene_list_tag} TopGenes on {gland_tag} Data", pad=10)
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)

    fname = f"{gene_list_tag}_TopGenes_on_{gland_tag}.pdf"
    g.savefig(fname, dpi=DPI, bbox_inches="tight")
    plt.close()
    return fname


# =============================== MAIN (已修改) ===============================
def main():
    # (修改) 1. 加载、拆分并获取细胞类型列表
    ads_dict, all_celltypes = load_and_split_anndata()
    if not ads_dict:
        print("[ERROR] 未能加载或拆分 AnnData。退出。")
        return
        
    top_dict = {}

    # (修改) 2. 按 Gland 筛选基因
    for gland in GLANDS:
        if gland not in ads_dict:
            continue
        print(f"\n[ {gland} ] 筛选腺体特异基因 …")
        
        # (修改) 传入 all_celltypes
        ct_deg = filter_gland_specific(gland, ads_dict, all_celltypes) 
        csvfile = save_deg_csv(ct_deg, gland)
        if csvfile:
            print(f"  CSV → {csvfile}")
            
        # (修改) 传入 all_celltypes
        top_dict[gland] = top_genes(ct_deg, all_celltypes) 
        print(f"  Top genes ({len(top_dict[gland])}) collected.")

    # (修改) 3. 绘图 (N x N 张)
    print("\n[ PLOTTING ] 开始绘制所有热图...")
    plot_count = 0
    for gene_gland, genes in top_dict.items(): 
        if not genes:
            print(f"  [INFO] {gene_gland} 基因列表为空，跳过绘图。")
            continue
        
        for data_gland, ad in ads_dict.items(): 
            
            # (修改) 传入 all_celltypes
            fname = plot_heatmap(genes, ad, data_gland, gene_gland, all_celltypes)
            if fname:
                print(f"  Heatmap saved: {fname}")
                plot_count += 1
            else:
                print(f"  Heatmap skipped: {gene_gland}_genes on {data_gland}_data")

    print(f"\n✅ 任务完成！共生成 {len(top_dict)} 个 CSV、{plot_count} 张 heatmap。")

if __name__ == "__main__":
    main()
