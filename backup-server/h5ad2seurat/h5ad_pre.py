import scanpy as  sc
import anndata as ad
#ea = ad.read_h5ad('./processed/allsc_sg_v9_scrublet.h5ad')
#ea = ad.read_h5ad('./original/allsc_sg_v9.h5ad')
ea = ad.read_h5ad('./processed/v9_all_log1p.h5ad')

if ea.raw is not None and ea.raw.X is not None:
    print('have raw data')
    adata = ad.AnnData(X=ea.raw.X.copy(),
                       var=ea.raw.var.copy(),
                       obs=ea.obs.copy())
else:
    print('no raw data')
sc.settings.figdir = ''
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
sc.settings.figdir = ''
sc.pl.umap(
        ea,
        color=["anno", "stage_new",'stage_old','leiden','gland','sample', "Pdgfra", "Acta2", "Col1a1",  # Mesenchymal
    "Ptprc","Cd8a", "Itgam",    # Immune
    "Fn1", "Mmp9", "Timp1" ],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save='v9_celltype.png',
    color_map='viridis'
    )


