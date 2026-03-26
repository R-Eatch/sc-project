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
sc.set_figure_params(
    dpi=100,          # 笔记本/预览分辨率
    dpi_save=300,     # 导出分辨率，期刊常用 300–600
    vector_friendly=True,  # 让散点在 PDF/SVG 中栅格化（默认就是 True）
    format='pdf'
)

# ========== Global Variables ==========

PCs = 10
res = 0.5
top_genes = 3000   # number of top variable genes
dataset = 'H-MG'  # sample name, will affect input/output file names
h5adpath = f'{dataset}_cleaned.h5ad'  # path to read the h5ad file

# Genes for feature-plot 
Featuregenes = [    "LPL",
    "CLDN3",
    "BTN1A1",
    "S100A1",
    "SPP1",
    "XBP1",
    "CLDN8"]

# Whether to plot feature genes
doFeatureplot = True

random.seed(2024)
np.random.seed(2024)


# In[1]:


def do_umap_plots(adata, sample_prefix="sample", feature_genes=None, do_feature=False):
    """
    Plot UMAP with basic annotations and optional feature genes.
    """
    sc.pl.umap(
        adata,
        color=["author_cell_type", "cell_type"],
        wspace=0.5,
        size=3,
        ncols=3,
        save=f'{sample_prefix}_basic_umap.png'
    )
    if do_feature and feature_genes is not None:
        
        sc.pl.umap(
            adata,
            color=feature_genes + ["author_cell_type", "cell_type"],
            color_map='viridis',
            ncols=2,
            frameon=False,
            legend_loc='on data',
            save=f"_{dataset}_featureplot.pdf",
            show=False,
            gene_symbols='feature_name'
        )
        print("Feature gene UMAP plots completed.")
        adata.var_names = adata.var["feature_name"]
        sc.pl.violin(
        adata,
        keys= feature_genes,
        groupby='author_cell_type', 
        jitter=0.4,
        rotation=90,
        size=2,stripplot=False,save=f"{sample_prefix}_celltype_violin.pdf",use_raw=False
        )
        sc.pl.stacked_violin(adata,var_names=feature_genes,gene_symbols="feature_name", groupby='author_cell_type', dendrogram=False,save=f"{sample_prefix}_violin_all.pdf")
    else:
        print("Skipping feature gene UMAP plots.")


# In[ ]:


# ========== Main Execution ==========

if __name__ == "__main__":
    print(f"Reading h5ad file: {h5adpath}")
    adata = sc.read_h5ad(h5adpath)
    do_umap_plots(
            adata=adata,
            sample_prefix=dataset,
            feature_genes=Featuregenes,
            do_feature=doFeatureplot
    )
# In[ ]:




