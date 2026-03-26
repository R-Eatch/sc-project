import pandas as pd
import scanpy as sc

# 修改为你的路径
pkl_in = "./v9_index_filter.pkl"
adata_in = "/data01/sunxuebo/project/scrnaseq/h5ad2seurat/processed/v9_all_counts.h5ad"

# 读取
df = pd.read_pickle(pkl_in)
adata = sc.read_h5ad(adata_in)

# 统一转为字符串比较（避免同值但 dtype 不同导致假不匹配）
ad_idx = pd.Index(adata.obs_names.astype(str), name="index")
pk_idx = pd.Index(df.index.astype(str), name="index")

# 统计
n_ad = ad_idx.size
n_pk = pk_idx.size
matched = ad_idx.intersection(pk_idx).size
unmatched_ad = ad_idx.difference(pk_idx).size
extra_pk = pk_idx.difference(ad_idx).size
dup_ad = n_ad - ad_idx.nunique()
dup_pk = n_pk - pk_idx.nunique()
same_order = (n_ad == n_pk) and ad_idx.equals(pk_idx)  # 是否长度一致且顺序完全相同

print({
    "adata_total": n_ad,
    "pkl_total": n_pk,
    "matched": matched,
    "unmatched_in_adata": unmatched_ad,
    "extra_in_pkl": extra_pk,
    "duplicate_index_in_adata": dup_ad,
    "duplicate_index_in_pkl": dup_pk,
    "coverage_%": round(matched / n_ad * 100, 2) if n_ad else 0.0,
    "same_order": same_order
})

# 如需查看少量未匹配项（仅调试用），解开下面两段注释：
# if unmatched_ad:
#     print("Unmatched (in adata only) head10:", list(ad_idx.difference(pk_idx))[:10])
# if extra_pk:
#     print("Extra (in PKL only) head10:", list(pk_idx.difference(ad_idx))[:10])

