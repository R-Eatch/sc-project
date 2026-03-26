#!/usr/bin/env python
# coding: utf-8

# In[110]:


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


# In[111]:


###global variable###


# In[112]:


PCs = 20
res = 0.3
top_genes=2000
dataset = 'S-AG'
h5adpath = f'../0.loom/{dataset}_for_subset.h5ad'
celltypelist = ['Luminal epithelial cells','Basal epithelial cells','Proliferating epithelial cells','Keratinized epithelial cells']
vars_use=['sample']


update_cell_type = True
clusterlist = ['']
#Featuregenes = ['Esr1','Epcam','Lalba','Top2a','Pgr','Prlr','Acta2','Elf5','Tcf4','Krt1','Ar','Pigr','Cd69','Adipoq','Lum','Vim','Ptprc','Lef1','Tpm2','Krt23','Krt10','Faah','Tm4sf1','Ppl','Wnt11','Krtdap','Sbsn','Dsp','Rab25','Aqp3','Shh','Atp1a1','Atp1b1','Procr']
random.seed(2024)
np.random.seed(2024)
cells = 5000
doFeatureplot = True
runharmony = True
subset_celltype = True
do_subset_cluster = False
random_subset_cells = False
use_scvi=False


###################function###############



def do_DEG(ea,dataset):
    #### Different gene test###
    sc.tl.rank_genes_groups(ea,groupby='leiden',method = 'wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ea,groupby='leiden',n_genes=5,save=f'{dataset}_dotplot_leiden-DEG.png',min_logfoldchange=0.25)
      # 提取 DEG 结果（每个 cluster 的 top genes）
    deg_df = sc.get.rank_genes_groups_df(ea, group=None,pval_cutoff=0.01,log2fc_min=0.25)  # group=None 返回所有组别拼接的 DataFrame
    deg_df.to_csv(f'{dataset}_ranked_genes_leiden-DEG.csv', index=False)


ea=sc.read_h5ad(f'{dataset}_for_DEG.h5ad')


do_DEG(ea=ea,dataset=dataset)





