#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

import os
# Disable hash randomization for reproducibility
os.environ['PYTHONHASHSEED'] = '0'

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import harmonypy as hm

print("Modules imported successfully.")
sc.settings.figdir = ''
# ========== Global Variables ==========

PCs = 10
res = 0.3
top_genes = 2000   # number of top variable genes
dataset = 'H-MG'  # sample name, will affect input/output file names
h5adpath = f'{dataset}.h5ad'  # path to read the h5ad file
vars_use = ['sample']  # variables used in Harmony

# List of clusters to keep (example). Adjust these to the cluster IDs you want to retain.
clusterlist = [2,4,5,11]

# Whether to update new cell type annotation
update_cell_type = False
NewCellType = {
    "StemCells": [2, 4, 7, 6, 10]
}

# Example: only keep certain original cell types
celltypelist = ['Luminal epithelial cells',
                'Basal epithelial cells',
                'Proliferating epithelial cells']
# Whether to subset by celltype or not
subset_celltype = False

# Genes for feature-plot (as in your original script)
Featuregenes = [
    'Esr1','Epcam','Top2a','Acta2','Elf5','Tcf4','Krt1','Prlr','Ar','Pigr',
    'Cd69','Adipoq','Lum','Pgr','Vim','Ptprc','Lalba','Lef1','Tpm2','Krt23',
    'Krt10','Faah','Tm4sf1','Ppl','Wnt11','Krtdap','Sbsn','Dsp','Rab25',
    'Aqp3','Shh','Atp1a1','Atp1b1','Procr','Wnt6','Krt5','Trp63','Lgr5','Krt5','Pip','Sox9','Sox10','Areg'
]

# Whether to plot feature genes
doFeatureplot = False

# Whether to run Harmony
runharmony = False

random.seed(2024)
np.random.seed(2024)


# In[1]:


def run_preprocess(adata, n_top_genes=3000):
    """
    Standard Scanpy preprocessing: normalize_total -> log1p -> highly_variable_genes.
    """
    adata.layer['counts']=adata.raw.X
    adata.X=adata.layer['counts']
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw = adata
    adata.layers["normalized"] = adata.X.copy()
    print("Data normalization and log1p completed.")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pl.highly_variable_genes(adata)
    print("Highly variable gene selection completed.")


def run_reduce_dimension(adata, run_harmony_flag=False, PCs=10, resolution=0.3, vars_use=None):
    """
    Scale data, PCA, neighbors, Leiden clustering, UMAP.
    Optionally run Harmony if run_harmony_flag is True.
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, log=False)
    print("PCA completed.")

    if run_harmony_flag and vars_use is not None:
        adata = run_harmony(adata, vars_use=vars_use)
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=PCs, random_state=42, n_neighbors=15)
        print("Neighbors graph constructed using Harmony output.")
    else:
        sc.pp.neighbors(adata, n_pcs=PCs, random_state=42)
        print("Neighbors graph constructed without Harmony.")

    sc.tl.leiden(adata, resolution=resolution, random_state=42)
    print("Leiden clustering done.")
    sc.tl.umap(adata, random_state=42)
    print("UMAP embedding finished.")

    return adata
def do_umap_plots(adata, sample_prefix="sample", feature_genes=None, do_feature=False):
    """
    Plot UMAP with basic annotations and optional feature genes.
    """
    sc.pl.umap(
        adata,
        color=["author_cell_type", "leiden"],
        wspace=0.5,
        size=3,
        ncols=3,
        save=f'{sample_prefix}_basic_umap.png'
    )
    if do_feature and feature_genes is not None:
        sc.pl.umap(
            adata,
            color=[*feature_genes,'author_cell_type','leiden'],
            wspace=0.5,
            size=3,
            ncols=4,
            save=f'{sample_prefix}_feature_umap.png'
        )
        print("Feature gene UMAP plots completed.")
    else:
        print("Skipping feature gene UMAP plots.")


def export_deg_result(adata):
    """
    Extract DE results from adata.uns['rank_genes_groups'] into a DataFrame.
    Returns only pval < 0.01 subset.
    """
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    data = []
    for group in groups:
        genes = result['names'][group]
        logfoldchanges = result['logfoldchanges'][group]
        pvals = result['pvals'][group]
        pvals_adj = result['pvals_adj'][group]
        for gene, lfc, pval, pval_adj in zip(genes, logfoldchanges, pvals, pvals_adj):
            data.append([group, gene, lfc, pval, pval_adj])
    df = pd.DataFrame(data, columns=['group', 'gene', 'logfoldchange', 'pval', 'pval_adj'])
    df1 = df[df['pval'] < 0.01]
    return df1


def do_deg(adata, sample_prefix="sample"):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        n_genes=5,
        save=f'{sample_prefix}_dotplot_leiden.png',
        min_logfoldchange=0.25
    )
    df1 = export_deg_result(adata)
    df1.to_csv(f'{sample_prefix}_ranked_genes_leiden.csv', index=False)
    print("DE analysis completed and results saved.")


# In[ ]:


# ========== Main Execution ==========

if __name__ == "__main__":
    print(f"Reading h5ad file: {h5adpath}")
    adata = sc.read_h5ad(h5adpath)
    # Preprocess
    run_preprocess(adata, n_top_genes=top_genes)
    adata1=adata[:,adata.var['highly_variable']]
    adata1.var['feature_name'].tolist()
    df1=pd.DataFrame(adata1.var['feature_name'].tolist(),columns=['gene'])
    df1.to_csv(f'{dataset}-HVGs.csv',index=False)
    # Run dimension reduction
    adata = run_reduce_dimension(
        adata,
        run_harmony_flag=runharmony,
        PCs=PCs,
        resolution=res,
        vars_use=vars_use
    )
    # UMAP plots
    do_umap_plots(
        adata=adata,
        sample_prefix=dataset,
        feature_genes=Featuregenes,
        do_feature=doFeatureplot
    )
    # DEG analysis
    do_deg(adata, sample_prefix=dataset)
    output_filename = f'{dataset}_cleaned.h5ad'
    adata.write(output_filename)
    print(f"All steps completed. The new h5ad file is saved as: {output_filename}")


# In[ ]:




