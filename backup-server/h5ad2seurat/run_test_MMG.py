import scanpy as  sc
import anndata as ad

adata = sc.read_h5ad('allcell_3species_v7.h5ad')

adata_subset=adata[(adata.obs['species'] == "M") & (adata.obs['gland'] == "MG")]

adata_subset.write('M-MG-3000.h5ad')
