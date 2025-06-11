#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as  sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import os
import tempfile
import scvi
import seaborn as sns
import torch
from rich import print
from scib_metrics.benchmark import Benchmarker


# In[ ]:


def run_scvi(adata,batch,res):
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    model.train()
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    return(adata)


# In[ ]:


###global variable###
PCs = 10 
res = 0.2
dataset1='R-MG'
dataset2='R-AG'
dataset = 'R-MAG'
path1=f"../../{dataset1}/1.subset/{dataset1}_cleaned.h5ad"
path2=f"../../{dataset2}/1.subset/{dataset2}_cleaned.h5ad"
vars_use = ['gland', 'stage']
random.seed(12345)
np.random.seed(12345)
NewCellType = {
    "StemCells": [0,5,6]
}
Featuregenes = ['Esr1','Epcam','Top2a','Acta2','Prlr','Tcf4','Ar','Lalba','Elf5','Gata3']
doFeatureplot = True
subset_celltype = True
re_find_hvgenes = True
random_state=2024
celltypelist = ['Lum-Immun-like','Lum-Ker-like','Lum-Fibro-like']
sc.settings.figdir = ''
#########################


# In[ ]:


adata1 = sc.read_h5ad(path1)
adata2 = sc.read_h5ad(path2)


# In[ ]:


adata1.obs = adata1.obs.drop(columns=["cellid"])
adata2.obs = adata2.obs.drop(columns=["cellid"])


# In[ ]:


if subset_celltype:
    adata1 = subsetcelltype(adata1,celltypelist)
    adata2 = subsetcelltype(adata2,celltypelist)
    unique_gland = '-'.join(adata1.obs['gland'].unique().tolist())
    adata1.obs['newcelltype'] = [unique_gland + '-' + str(newtype) for newtype in adata1.obs['newcelltype']]
    unique_gland = '-'.join(adata2.obs['gland'].unique().tolist())
    adata2.obs['newcelltype'] = [unique_gland + '-' + str(newtype) for newtype in adata2.obs['newcelltype']]


# In[ ]:


adata_combined = adata1.concatenate(adata2,batch_key = 'glandbatch')


# In[ ]:


if re_find_hvgenes:
    adata_combined.X = adata_combined.layers['normalized']
    # Normalizing to median total counts
    sc.pp.normalize_total(adata_combined,target_sum= 1e4)
    # Logarithmize the data
    sc.pp.log1p(adata_combined)
    print('Finish normalized')
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000)
    print('Finish Varible genes')
    #sc.pp.scale(adata_combined,max_value=10)


# In[ ]:


run_scvi(adata=adata_combined,batch='glandbatch',res=res)


# In[ ]:


sc.pp.neighbors(adata_combined,use_rep=SCVI_LATENT_KEY,random_state=random_state)
sc.tl.leiden(adata_combined,resolution = res,random_state=random_state)
print('Finish clustering')
sc.tl.umap(adata_combined,random_state=random_state)
print('Finish UMAP')


# In[ ]:


sc.pl.umap(
        adata_combined,
        color=["newcelltype", "celltype", "stage",'sample','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_merged_celltype.png',
    color_map='viridis'
    )


# In[ ]:


adata_combined.write(f'{dataset}_cleaned.h5ad')

