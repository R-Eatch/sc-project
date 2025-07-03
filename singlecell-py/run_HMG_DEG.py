#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as  sc
import anndata as ad
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
print('import successful')


# In[ ]:


sc.settings.figdir = ''
# ========== Global Variables ==========
dataset = 'H-MG'  # sample name, will affect input/output file names
h5adpath = f'../{dataset}_cleaned.h5ad'  # path to read the h5ad file


# In[ ]:


adata=sc.read_h5ad(h5adpath)


# In[7]:


adata


# In[5]:


def do_deg(adata, sample_prefix="sample",group='leiden'):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby=group, method='wilcoxon',pts=True)
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=group,
        n_genes=5,
        save=f'{sample_prefix}_dotplot_{group}.png',
        min_logfoldchange=0.25
    )
    df1 = sc.get.rank_genes_groups_df(adata, gene_symbols='feature_name')
    df1.to_csv(f'{sample_prefix}_ranked_genes_{group}.csv', index=False)
    print("DE analysis completed and results saved.")


# In[ ]:


do_deg(adata,sample_prefix=dataset,group='author_cell_type')


# In[ ]:


do_deg(adata,sample_prefix=dataset,group='cell_type')


# In[ ]:




