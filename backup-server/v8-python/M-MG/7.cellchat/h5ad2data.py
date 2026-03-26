import scanpy as  sc
import anndata as ad

dataset= 'M-MG'
ea = ad.read_h5ad(f'../1.subset/{dataset}_cleaned.h5ad')

ea.X=ea.layers['normalized']
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
ea.write(f'{dataset}_raw.h5ad')

