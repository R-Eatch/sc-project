#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')

import numpy as np

from matplotlib import pyplot as plt

from pyslingshot import Slingshot
import scanpy as sc


# In[34]:


start_node=''
end_nodes=[]
celltype_key='newcelltype'
obsm_key='X_umap'
dataset="M-MG"


# In[3]:


#adata=sc.read_h5ad(f"D:/111/{dataset}_cleaned.h5ad")


# In[3]:


adata=sc.read_h5ad(f"../1.subset/{dataset}_cleaned.h5ad")


# In[4]:


adata


# In[30]:


if start_node and end_nodes:
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose',end_nodes=end_nodes,start_node=start_node)
elif start_node:
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose',start_node=start_node)
else:
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose')
    


# In[35]:


fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 10))
custom_xlim = (-12, 12)
custom_ylim = (-12, 12)
# plt.setp(axes, xlim=custom_xlim, ylim=custom_ylim)

if start_node and end_nodes:
    print(f'start:{start_node} and end:{end_nodes}')
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose',end_nodes=end_nodes,start_node=start_node)
elif start_node:
    print(f'start:{start_node}')
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose',start_node=start_node)
else:
    print('666')
    slingshot = Slingshot(adata, celltype_key=celltype_key, obsm_key=obsm_key, debug_level='verbose')
slingshot.fit(num_epochs=1, debug_axes=axes)
fig.savefig(f'{dataset}_slingshot1.png',dpi=300)
fig.savefig(f'{dataset}_slingshot1.pdf',dpi=300)
fig.show()


# In[41]:


fig, axes = plt.subplots(ncols=2, figsize=(12, 4))
axes[0].set_title('cluster')
axes[1].set_title('Pseudotime')
slingshot.plotter.curves(axes[0], slingshot.curves)
slingshot.plotter.clusters(axes[0], labels=np.arange(slingshot.num_clusters), s=4, alpha=0.5)
slingshot.plotter.curves(axes[1], slingshot.curves)
slingshot.plotter.clusters(axes[1], color_mode='pseudotime', s=5)

axes[0].legend(loc='upper left', bbox_to_anchor=(1.05, 1))
axes[1].legend(loc='upper left', bbox_to_anchor=(1.05, 1))

# 自动调整子图排版，防止图例遮挡
fig.tight_layout()
fig.savefig(f'{dataset}_slingshot2.png',dpi=300)
fig.savefig(f'{dataset}_slingshot2.pdf',dpi=300)
fig.show()


# In[10]:


# NOTE: the Slingshot class has a property which has the pseudotime that is used to 
# color the plot above
pseudotime = slingshot.unified_pseudotime


# In[11]:


pseudotime


# In[ ]:




