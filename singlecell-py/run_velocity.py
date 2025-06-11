#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import plotly.express as px
import matplotlib.pyplot as plt
print('import successful')


# In[2]:


import scvelo as scv
scv.logging.print_version()
scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[3]:


###global variable###
dataset = 'M-MG'
ad_path = f'../1.subset/{dataset}_cleaned.h5ad'
use_default_method=True


# In[4]:


###funciton###


# In[5]:


def run_plotly(adata):

    umap_data = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
    umap_data['cell_ids'] = adata.obs.index 
    umap_data['time'] = adata.obs['time'] 

    fig = px.scatter(umap_data, x='UMAP1', y='UMAP2', hover_data=['cell_ids', 'time'],width=800,height=450)
    fig.show()
    fig.write_html('interactive_umap_plot.html')


# In[4]:


ad_path='D:/111/M-MG_for_subset.h5ad'


# In[7]:


adata1 = ad.read(ad_path)


# In[9]:


run_plotly(adata1)


# In[ ]:


adata = ad.read(ad_path)
print("Before normalization:")
print("X shape:", adata.X.shape)  # 查看基因表达矩阵的形状

scv.pp.filter_and_normalize(adata1,log=False, n_top_genes=None)
print("\nAfter normalization:")
print("X shape:", adata.X.shape)  # 基因表达矩阵不应变化


# In[10]:


scv.pp.filter_and_normalize(adata1,log=False, n_top_genes=None)
sc.pp.neighbors(adata1, n_pcs=10, n_neighbors=30)
scv.pp.moments(adata1, n_pcs=10, n_neighbors=30)


# In[11]:


### run velocity ###
if use_default_method:    
    scv.tl.velocity(adata1,mode='stochastic')
    print("use default mode =  stochastic ")
else:
    scv.tl.recover_dynamics(adata1)
    scv.tl.velocity(adata1, mode='dynamical')
    print("use mode =  dynamic ")


# In[13]:


scv.tl.velocity_graph(adata1)


# In[14]:


scv.pl.velocity_embedding_stream(adata1, basis='umap',color='newcelltype',save=f'./{dataset}_velo.png')


# In[15]:


scv.pl.proportions(adata1,save=f'./{dataset}_proportions.png')


# In[18]:


adata1.write(f'{dataset}_velofinished.h5ad')
print(f"save velofinished file : {dataset}_velofinished.h5ad")


# In[16]:


scv.tl.score_genes_cell_cycle(adata1)
scv.pl.scatter(adata1, 
               color_gradients=['S_score', 'G2M_score'], 
               palette=['green', 'orange'], 
               smooth=True, perc=[5, 90],
               save=f'./{dataset}_cycle.png'
               )


# In[17]:


scv.tl.velocity_pseudotime(adata1)
scv.pl.scatter(adata1, color='velocity_pseudotime', cmap='gnuplot',save=f'./{dataset}_velo_pseu.png')


# In[19]:


###test code do not run###


# In[ ]:




