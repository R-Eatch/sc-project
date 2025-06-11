#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
# Set the environment variable to disable hash randomization
os.environ['PYTHONHASHSEED'] = '0'
#If you find discrepancies between the UMAP and Leiden results, indicating reproducibility issues, 
#please refer to this link: https://github.com/scverse/scanpy/issues/1009 for potential solutions.
#https://umap-learn.readthedocs.io/en/latest/reproducibility.html
# set export OMP_NUM_THREADS=1 Multithreading can lead to reproducibility issues!!!
#Theoretically, minor differences in UMAP and inconsistencies in Leiden clustering do not affect their biological significance.
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


# In[ ]:


PCs = 10
res = 0.3
top_genes=2000
datasetlist=['M-MG','R-MG','S-MG','R-AG','R-CG','S-AG']
vars_use=['sample']
input_file_path = "/data01/sunxuebo/project/scrnaseq/h5ad2seurat/processed/v9_all_counts.h5ad"
Featuregenes = ['Esr1','Epcam','Lalba','Top2a','Pgr','Prlr','Acta2','Elf5','Tcf4','Krt1','Ar','Pigr','Cd69','Adipoq','Lum','Vim','Ptprc','Lef1','Tpm2','Krt23','Krt10','Faah','Tm4sf1','Ppl','Wnt11','Krtdap','Sbsn','Dsp','Rab25','Aqp3','Shh','Atp1a1','Atp1b1','Procr']
random.seed(2024)
np.random.seed(2024)
doFeatureplot = True
runharmony = True
use_scvi=False
celltypelist=['Fibroblasts','B lymphocytes','T lumphocytes','Innate immune cells']


# In[ ]:


def split_h5ad_file(input_file):
    # Read the input h5ad file
    adata = ad.read_h5ad(input_file)
    # 1) Create a new/overwrite species column from the first letter of 'sample'
    adata.obs['species'] = adata.obs['sample'].str[0]
    # 2) Define conditions for splitting based on species and gland
    #    (including the two new groups 'S-SG' and 'M-SG')
    conditions = {
        "M-MG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "MG"),
        "R-MG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "MG"),
        "R-AG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "AG"),
        "S-MG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "MG"),
        "S-AG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "AG"),
        "R-CG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "CG")
    }
    # Base directory for the output files
    base_dir = "/data01/sunxuebo/project/scrnaseq/v9-immune"

    # Process each condition and save the subset
    for group_name, condition in conditions.items():
        subset_adata = adata[condition].copy()
        group_path = os.path.join(base_dir, group_name)
        if not os.path.exists(group_path):
            os.makedirs(group_path)
        subset_adata.write(os.path.join(group_path, f"{group_name}-immune.h5ad"))
        print(f"Saved {group_name}.h5ad to {group_path}")


# In[ ]:


def run_harmony(adata,vars_use=vars_use):
    print('running harmony')
    pca_result = adata.obsm['X_pca']
    ho = hm.run_harmony(pca_result, adata.obs, vars_use,random_state=42)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print('finished harmony')
    return adata


# In[ ]:


def run_preprocess(adata,top_genes):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    adata.raw=adata
    adata.layers["normalized"] = adata.X.copy()
    print('Finish normalized')
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
    sc.pl.highly_variable_genes(adata)
    print('Finish Varible genes')


# In[ ]:


def run_reduceDimension(ea,use_scvi,runharmony,PCs,res):
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
        run_harmony(ea)
        sc.pp.neighbors(ea, use_rep='X_pca_harmony',n_pcs=PCs,random_state=42,n_neighbors=15)
        print('Finish neighbors')
    else:
        sc.pp.neighbors(ea,n_pcs=PCs,random_state=42)
        print('Finish neighbors')
    sc.tl.leiden(ea,resolution = res,random_state=42)
    print('Finish clustering')
    sc.tl.umap(ea,random_state=42)
    print('Finish UMAP')


# In[ ]:


def do_umap_plots(ea,dataset,Featuregenes,doFeatureplot):
    sc.settings.figdir = f'./{dataset}'
    sc.pl.umap(
        ea,
        color=["newcelltype", "anno", "stage",'sample','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_cluster.png',
    color_map='viridis'
    )

    if doFeatureplot:
        sc.pl.umap(
            ea,
            color=[*Featuregenes,'newcelltype',"stage",'sample','leiden','gland'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
        ncols=4,
        save=f'{dataset}_cluster_feature.png',
        color_map='plasma'
        )
        print('Featureplot Finished')
    else:
        print('666666')


# In[2]:


def update_adata(adata2,dataset):
    adata2.obs['newcelltype']=adata2.obs['anno']
    adata2.obs['celltype']=adata2.obs['anno']
    adata2.obs['stage']=adata2.obs['stage_new']
    adata2.layers['counts']=adata2.X.copy()
    adata2.obs['dataset']=dataset
    adata2.obs['cellid'] = (adata2.obs['sample'].astype(str) + '-' + adata2.obs_names.astype(str)).str.replace(r'(-\d+)+$', '', regex=True)
    return adata2


# In[ ]:


# Call the function with the direct input file path
split_h5ad_file(input_file_path)


# In[ ]:


for dataset in datasetlist:
    adata1_path=f'/data01/sunxuebo/project/scrnaseq/v8-python/{dataset}/1.subset/{dataset}_cleaned.h5ad'
    adata2_path=f'./{dataset}/{dataset}-immune.h5ad'
    adata1=sc.read_h5ad(adata1_path)
    adata2_raw=sc.read_h5ad(adata2_path)
    adata2_raw=update_adata(adata2=adata2_raw,dataset=dataset)
    adata2=adata2_raw[adata2_raw.obs['newcelltype'].isin(celltypelist)]
    adata = adata1.concatenate(adata2,batch_key = 'cellbatch')
    run_preprocess(adata,top_genes)
    run_reduceDimension(adata,use_scvi,runharmony,PCs,res)
    do_umap_plots(adata=adata,Featuregenes=Featuregenes,doFeatureplot=doFeatureplot)
    adata.write(f'./{dataset}/{dataset}_lum-imm.h5ad')
    adata.X=adata.layers['normalized']
    adata.write(f'./{dataset}/{dataset}_lum-imm_raw.h5ad') ## for h5ad2seurat and cellchat2
    del adata1
    del adata2
    del adata2_raw
    del adata

