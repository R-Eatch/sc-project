#!/usr/bin/env python
# coding: utf-8

# # Erythroid lineage

# ## 1. Loading packages

# In[2]:


import sys
sys.path = [".."] + sys.path # cd to the velovgi project path

import numpy as np
import matplotlib.pyplot as plt

import scvelo as scv
import velovgi

from pytorch_lightning import loggers
from torch_geometric import seed_everything

seed = 0
seed_everything(seed)


# ## 2. Loading anndata object

# The ideal estimated direction of differentiation: Blood progenitors 1 -> Blood progenitors 2 -> Erythroid 1 -> Erythroid 2 -> Erythroid 3 

# In[21]:


adata = scv.datasets.gastrulation_erythroid("../1.subset/S-MAG_cleaned.h5ad")
batch_key = "stage" # the column name identifing batch in obs
cluster_key = "newcelltype" # the column name identifing cluster in obs
scv.pl.umap(adata, color=[cluster_key, batch_key], legend_loc="right")
adata


# ## 3. Preprocessing Data

# The preprocessing process includes:
# - quality control
# - normalization
# - hvg(highly variable gene) selection
# - pca
# - neighborhood construction between adjacent batches in time series
# - moment

# In[5]:


subsample_adata = velovgi.pp.preprocess(adata, sample_mode="random", batch_key=batch_key)


# ## 4. Training VeloVGI model and export result

# VeloVGI can be trained with scvi-tools styles.

# In[6]:


logger = loggers.TensorBoardLogger(save_dir="./log", name="base") # TensorBoard logging file
velovgi.tl.VELOVGI.setup_anndata(adata=subsample_adata, spliced_layer="Ms", unspliced_layer="Mu")
velovgi_model = velovgi.tl.VELOVGI(subsample_adata)
velovgi_model.train(logger=logger)


# The  parameters related to RNA velocity can be exported from VeloVGI model to Anndata object, such as transcription rates $\alpha$, splicing rates $\beta$, degradation rates $\gamma$.

# In[7]:


velovgi.tl.add_velovi_outputs_to_adata(subsample_adata, velovgi_model) # 模型输出
velovgi.pp.moment_recover(adata, subsample_adata) # 恢复
subsample_adata, adata


# ## 5. Visualization

# ### 5.1 Stream plot

# RNA velocity stream plot can be visualized.

# In[8]:


scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, color=[cluster_key, batch_key], legend_loc="right",save=f'./{dataset}_veloVGI.png')


# ### 5.2 Hierarchical embedding

# Hierarchical embedding  shows heterogeneity of cell types in each batch.

# In[10]:


figsize=(10,10)
fig, ax = plt.subplots(figsize=figsize)

transform_matrix_params = dict(p=0.5, q=2 , theta=np.pi/20) # transition matrix parameters
sep = 25 # hierarchical gap parameters
embedding_plot_params = dict() # scatter plot parameters

velovgi.pl.draw_batch_layer_embedding(adata, cluster_key=cluster_key, batch_key=batch_key, transform_matrix_params=transform_matrix_params, sep=sep, embedding_plot_params=embedding_plot_params, ax=ax)


# The 3D interactive graph allow for more visualization of batch-to-batch differences.

# In[11]:


velovgi.pl.draw_batch_layer_embedding_3d(adata, cluster_key=cluster_key, batch_key=batch_key, size=3)


# ### 5.3 Neighborhood visualization

# The batches are arranged clockwise on the circos plot, and the thickness of the bands between batches represents the number of neighbors between batches.

# In[13]:


velovgi.pl.draw_batch_circos_ax(adata, batch_key=batch_key)


# The number of internal (knn) neighbors per batch and the number of external (bnn) neighbors between batches are shown numerically in the umap plot, respectively.

# In[14]:


velovgi.pl.draw_batch_nn_umap(adata, batch_key=batch_key, cluster_key=cluster_key)


# ## 6. Saving Result

# 1. Anndata saving.

# In[16]:


adata_dir = "./data/adata"
velovgi.tl.write_adata(adata, adata_dir)


# In[18]:


adata = velovgi.tl.read_adata(adata_dir)
adata


# 2. Model model saving. 

# In[19]:


model_dir = "./model/base"
velovgi_model.save(model_dir)

