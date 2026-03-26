#!/usr/bin/env python
# coding: utf-8

# In[2]:

from __future__ import annotations

from typing import List

import pandas as pd
from anndata import AnnData

__all__ = ["filter_stage_top_genes"]


def filter_stage_top_genes(
    adata: AnnData,
    csv_path: str,
    top_n: int = 100,
    *,
    group_col: str = "group",
    gene_col: str = "gene",
    inplace: bool = False,
) -> AnnData:
    """Return an AnnData containing only the top‑*n* genes per stage.

    Parameters
    ----------
    adata
        The input AnnData object.
    csv_path
        Path to the CSV file already sorted by significance.
    top_n
        Number of most‑significant genes to keep per stage (default 100).
    group_col, gene_col
        Column names in the CSV for stage label and gene symbol.
    inplace
        If *True*, modifies *adata* in‑place and returns it; otherwise
        returns a **copy** with only the selected genes.
    """
    # 1. Read minimal columns only
    df = pd.read_csv(csv_path, usecols=[group_col, gene_col])

    # 2. Grab top_n rows per stage and build the union gene set
    top_genes: List[str] = (
        df.groupby(group_col).head(top_n)[gene_col].unique().tolist()
    )

    # 3. Subset AnnData
    mask = adata.var_names.isin(top_genes)
    if inplace:
        adata._inplace_subset_var(mask)
        return adata
    return adata[:, mask].copy()

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix
import os


# In[3]:


#A=pd.read_csv('D:/111/maps/msrb/ms_to_rb.txt',sep='\t',index_col=0,header=None)
#A.head()


# In[ ]:


####function#####


# In[4]:


def cut_version_number(sp1,sp2):
    df=pd.read_csv(f'./maps/{sp1}{sp2}/{sp1}_to_{sp2}.txt',sep='\t',index_col=0,header=None)
    df.index=df.index.str.replace(r'\.\d+$','',regex=True)
    df[1]=df[1].str.replace(r'\.\d+$','',regex=True)
    df.to_csv(f'./maps/{sp1}{sp2}/{sp1}_to_{sp2}.txt',sep='\t',header=None)
    
    df=pd.read_csv(f'./maps/{sp1}{sp2}/{sp2}_to_{sp1}.txt',sep='\t',index_col=0,header=None)
    df.index=df.index.str.replace(r'\.\d+$','',regex=True)
    df[1]=df[1].str.replace(r'\.\d+$','',regex=True)
    df.to_csv(f'./maps/{sp1}{sp2}/{sp2}_to_{sp1}.txt',sep='\t',header=None)
    return None


# In[62]:


def dropobs(adata,obslist,sp):
    for i in adata.obs.columns:
        if i not in obslist:
            del adata.obs[i]
        else:
            adata.obs[i] = adata.obs[i].astype('category')
    del adata.var
    del adata.varm
    del adata.obsm
    del adata.uns
    del adata.obsp
    print(adata)
    return adata


# In[63]:


def clean_h5ad(dataset,df,sp,obslist,use_DEG):
    adata=sc.read_h5ad(f"/data01/sunxuebo/project/scrnaseq/v8-python/{dataset}/1.subset/{dataset}_cleaned.h5ad")
    if use_DEG:
        adata=filter_stage_top_genes(adata=adata,csv_path=DEG_PATH)
    adata.raw = None
    if sp == 'ms':
        mapping=dict(zip(df["Gene name"],df["Gene stable ID"]))
    if sp =='rb':
        mapping=dict(zip(df["Gene name"],df["Rabbit gene stable ID"]))
    if sp=='sg':
        mapping=dict(zip(df["Gene name"],df["sg"]))
    adata.X=adata.layers['counts']
    adata.X = csr_matrix(adata.X)
    adata.var_names = adata.var_names.map(mapping)
    del adata.layers
    adata=dropobs(adata,obslist,sp=sp)
    if use_DEG:
        adata.write(f"./{dataset}_counts_DEG.h5ad")
    else:
        adata.write(f"./{dataset}_counts.h5ad")
    print(adata)
    return None   


# In[ ]:


def preprocess_h5ad(dataset1, dataset2, df, sp1, sp2, obslist,use_DEG):
    if use_DEG:
        fn1 = f'./{dataset1}_counts_DEG.h5ad'
        fn2 = f'./{dataset2}_counts_DEG.h5ad'
    else:
        fn1 = f'./{dataset1}_counts.h5ad'
        fn2 = f'./{dataset2}_counts.h5ad'
    # 如果 fn1 和 fn2 都存在，跳过整个函数
    if os.path.exists(fn1) and os.path.exists(fn2):
        print(f"Skipping preprocessing: Both files already exist: {fn1}, {fn2}")
        return  # 直接返回，跳过后续的操作

    # 清理数据
    clean_h5ad(dataset=dataset1, df=df, sp=sp1, obslist=obslist,use_DEG=use_DEG)
    clean_h5ad(dataset=dataset2, df=df, sp=sp2, obslist=obslist,use_DEG=use_DEG)

    ##### SAMAP - generate SAM file #############
    resolutions = {sp1: 0.3, sp2: 0.3}
    filenames = {sp1: fn1, sp2: fn2}
    sams = {sp1: fn1, sp2: fn2}

    # 执行 SAMAP
    sm = SAMAP(
        sams,
        f_maps='/data01/sunxuebo/project/scrnaseq/v8-python/samap/re-maps/',
        save_processed=True,
        resolutions=resolutions
    )


# In[ ]:


###########global variable############


# In[23]:

use_DEG=True
gene_map_csv = '/data01/sunxuebo/project/scrnaseq/v8-python/data/ortho_mrs.csv'  
obslist=['newcelltype','celltype','gland','species','stage','sample','leiden']
cut_df=False ## first run set true
datasetlist=['M-MG','R-MG','S-MG','R-AG','R-CG','S-AG']
#datasetlist=['M-SG','S-SG']
processed_pairs=set()
DEG_PATH='/data01/sunxuebo/project/scrnaseq/v8-python/MRS-MG/3.harmony/ALL-MG_stage_genes.csv'

# In[3]:


df=pd.read_csv("/data01/sunxuebo/project/scrnaseq/v8-python/data/ortho_mrs.csv")
#df=pd.read_csv('D:/111/ortho_mrs.csv')
df.head()

# In[9]:


####adjust data format######


# In[11]:


if cut_df:
    cut_version_number(sp1,sp2)
    print("cut finfished")


# In[ ]:


for dataset1 in datasetlist:
    for dataset2 in datasetlist:  
        if dataset1:
            ls=dataset1.split(sep='-')
        if ls[0]=='M':
            sp1='ms'
        elif ls[0]=='R':
            sp1='rb'
        elif ls[0]=='S':
            sp1='sg'
        gd1=ls[1]
        if dataset2:
            ls=dataset2.split(sep='-')
        if ls[0]=='M':
            sp2='ms'
        elif ls[0]=='R':
            sp2='rb'
        elif ls[0]=='S':
            sp2='sg'
        gd2=ls[1] 
        if dataset1 == dataset2:
            continue
        pair = tuple(sorted([dataset1, dataset2]))
        if pair in processed_pairs:
            continue
        preprocess_h5ad(dataset1,dataset2,df,sp1,sp2,obslist,use_DEG)
        processed_pairs.add(pair)
        print(pair)

