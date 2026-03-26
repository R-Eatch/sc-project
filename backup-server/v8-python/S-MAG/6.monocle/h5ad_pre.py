import scanpy as  sc
import anndata as ad

dataset= 'S-MG'
ea = ad.read_h5ad(f'../1.subset/{dataset}_cleaned.h5ad')
if ea.raw is not None and ea.raw.X is not None:
    print('have raw data')
    adata = ad.AnnData(X=ea.raw.X.copy(),
                       var=ea.raw.var.copy(),
                       obs=ea.obs.copy())
else:
    print('no raw data')
# 查看表达矩阵 X 的形状
print("X shape:", ea.shape)

# 查看细胞元数据（obs）的前5行
print("\n=== Observations (obs) head ===")
print(ea.obs.head())

# 查看基因元数据（var）的前5行
print("\n=== Variables (var) head ===")
print(ea.var.head())

# 查看 AnnData 对象的整体结构
print("\n=== AnnData object summary ===")
print(ea)
adata.write(f'{dataset}_raw.h5ad')

