#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
os.environ['PYTHONHASHSEED'] = '0'
os.environ['OPENBLAS_NUM_THREADS'] = '8'
os.environ['OMP_NUM_THREADS'] = '8'
os.environ['MKL_NUM_THREADS'] = '8'
import scanpy as  sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import loompy
import harmonypy as hm
print('import successful')


# In[2]:


PCs = 50
res = 0.2
top_genes=3000
glandlist=['MG','SG','EG']
vars_use=['sample']
Featuregenes = ['Esr1','Epcam','Lalba','Top2a','Pgr','Prlr','Acta2','Elf5','Tcf4','Krt1','Ar','Pigr','Cd69','Adipoq','Lum','Vim','Ptprc','Lef1','Tpm2','Krt23','Krt10','Faah','Tm4sf1','Ppl','Wnt11','Krtdap','Sbsn','Dsp','Rab25','Aqp3','Shh','Atp1a1','Atp1b1','Procr','Col1a1','Krt8','Krt18','Krt14','Krt7','Krt5','Csn1s1','Csn2','Csn3','Xdh']
random.seed(2024)
np.random.seed(2024)
doFeatureplot = False
runharmony = True
use_scvi=False
###0313
PCs = 20
res = 0.2
top_genes=2000
glandlist=['MG','SG','EG']
vars_use=['sample']
###0313

# In[3]:


def run_harmony(adata,vars_use=vars_use):
    print('running harmony')
    pca_result = adata.obsm['X_pca']
    ho = hm.run_harmony(pca_result, adata.obs, vars_use,random_state=42)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print('finished harmony')
    return adata


# In[4]:


def run_preprocess(adata,top_genes):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    adata.raw = adata
    adata.layers["normalized"] = adata.X.copy()
    print('Finish normalized')
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
    sc.pl.highly_variable_genes(adata)
    print('Finish Varible genes')


# In[5]:


def run_reduceDimension(ea,use_scvi,runharmony,PCs,res,vars_use):
    import scanpy as  sc
    sc.pp.scale(ea,max_value=10)
    sc.tl.pca(ea,mask_var="highly_variable")
    sc.pl.pca_variance_ratio(ea,log = False)
    if use_scvi:   
        import os
        import tempfile
        import scanpy as sc
        import scvi
        import seaborn as sns
        import torch
        from scib_metrics.benchmark import Benchmarker
        scvi.settings.seed = 0
        print("Last run with scvi-tools version:", scvi.__version__)
        ea=run_scvi(ea)
        sc.pl.embedding(
            ea,
            basis=SCVI_MDE_KEY,
            color=['stage', "leiden",'celltype'],
            frameon=False,
            ncols=1,
            save=f'{dataset}_SCVI.png',
        )
    else:
        print('skip scvi')
    if runharmony:
        run_harmony(ea,vars_use)
        sc.pp.neighbors(ea, use_rep='X_pca_harmony',n_pcs=PCs,random_state=42,n_neighbors=15)
        print('Finish neighbors')
    else:
        sc.pp.neighbors(ea,n_pcs=PCs,random_state=42)
        print('Finish neighbors')
    sc.tl.leiden(ea,resolution = res,random_state=42)
    print('Finish clustering')
    #ea.obsm['X_umap_original'] = ea.obsm['X_umap'].copy()
    sc.tl.umap(ea,random_state=42)
    ea.obsm['X_umap_original'] = ea.obsm['X_umap'].copy()
    print('Finish UMAP')


# In[6]:


def do_umap_plots(ea,dataset,Featuregenes,doFeatureplot):
    #os.makedirs('./figure', exist_ok=True)
    #sc.settings.figdir = './figure'
    os.makedirs('./figure-0313', exist_ok=True)
    sc.settings.figdir = './figure-0313'

    sc.pl.umap(
        ea,
        color=["celltype",'leiden','stage'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_cluster.png',
    color_map='viridis',
    )
    sc.pl.embedding(
        ea, basis='umap_original',
        color=["celltype",'leiden','stage'],
        wspace=0.5,
        size=3,
        ncols=3,
        save=f'{dataset}_cluster_ori.png',
        color_map='viridis')
    if doFeatureplot:
        sc.pl.umap(
            ea,
            color=[*Featuregenes,'celltype'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
            ncols=4,
            save=f'{dataset}_cluster_feature.png',#gene_symbols='gene_name',
            color_map='plasma'
        )
        sc.pl.embedding(
            ea,basis='umap_original',
            color=[*Featuregenes,'celltype'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
            ncols=4,
            save=f'{dataset}_cluster_feature_ori.png',#gene_symbols='gene_name',
            color_map='viridis')
        print('Featureplot Finished')
    else:
        print('666666')


# In[ ]:


def main():
    for dataset in glandlist:
        print(f'processing {dataset}')
        adata=sc.read_h5ad(f'./data/Final_{dataset}_merged.h5ad')
        run_preprocess(adata,top_genes=top_genes)
        run_reduceDimension(adata,use_scvi,runharmony,PCs,res,vars_use)
        do_umap_plots(adata,dataset,Featuregenes,doFeatureplot)
        #adata.write(f'./data/{dataset}_pp.h5ad')
        adata.write(f'./data/{dataset}_pp-0313.h5ad')
        print(f'{dataset} finished')


# In[ ]:


if __name__ == "__main__":
    main()

