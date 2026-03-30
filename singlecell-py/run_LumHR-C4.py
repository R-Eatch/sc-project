#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as  sc
import anndata as ad
import loompy
import scvelo as scv
import scanpy as sc
import os
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.collections as mc
print('import successful')


# In[2]:


dataset = 'S-MG'


# In[ ]:


adata=sc.read_h5ad(f"../../S-MG/1.subset/{dataset}_cleaned.h5ad")
adata


# In[4]:


adata_sub = adata[adata.obs['newcelltype'].isin(['LumHR'])].copy()
sc.tl.leiden(adata_sub,resolution=0.3)
sc.pl.umap(adata_sub,color='leiden',show=False)


# In[5]:


fig = plt.gcf()
for ax in fig.axes:
    for coll in ax.collections:
        if isinstance(coll, mc.PathCollection):   # 只栅格化散点
            coll.set_rasterized(True)
# 3) 保存为 PDF：dpi 只影响被 rasterize 的那部分（也就是点）
fig.savefig(f'{dataset}-umap-lumhr-c4.pdf', bbox_inches="tight", dpi=300)
plt.close(fig)


# In[6]:


sc.tl.rank_genes_groups(adata_sub,groupby='leiden',key_added='subHR')
sc.tl.dendrogram(adata_sub, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(adata_sub,groupby='leiden',n_genes=10,save=f'{dataset}-LumHR_dotplot_leiden.pdf',min_logfoldchange=0.25,key='subHR')


# In[7]:


df=sc.get.rank_genes_groups_df(adata_sub,key='subHR',group=None,log2fc_min=0.25,pval_cutoff=0.01)
df.to_csv('S-MG-LumHR.csv',index=False)
df1=df.groupby('group').head(50)
glist=df1.loc[df1['group']=='4','names'].tolist()
exclude_prefix = ('Rpl', 'Rps')
filtered_list = [gene for gene in glist if not gene.startswith(exclude_prefix)]
filtered_list


# In[8]:


sc.tl.score_genes(adata,gene_list=filtered_list,score_name='c4')


# In[9]:


sc.pl.umap(adata,color=[*filtered_list,'c4','leiden'],ncols=3,show=False,cmap="coolwarm")
fig = plt.gcf()
for ax in fig.axes:
    for coll in ax.collections:
        if isinstance(coll, mc.PathCollection):   # 只栅格化散点
            coll.set_rasterized(True)
# 3) 保存为 PDF：dpi 只影响被 rasterize 的那部分（也就是点）
fig.savefig(f'{dataset}-umap-lumhr-c4.pdf', bbox_inches="tight", dpi=300)
plt.close(fig)


# In[ ]:





# In[ ]:


dataset = 'S-AG'


# In[ ]:


adata2=sc.read_h5ad(f"../../S-AG/1.subset/{dataset}_cleaned.h5ad")


# In[ ]:


sc.tl.score_genes(adata2,gene_list=filtered_list,score_name='c4')


# In[ ]:


sc.pl.umap(adata2,color=[*filtered_list,'c4','leiden'],ncols=3,show=False,cmap="coolwarm")
fig = plt.gcf()
for ax in fig.axes:
    for coll in ax.collections:
        if isinstance(coll, mc.PathCollection):   # 只栅格化散点
            coll.set_rasterized(True)
# 3) 保存为 PDF：dpi 只影响被 rasterize 的那部分（也就是点）
fig.savefig(f'{dataset}-umap-lumhr-c4', bbox_inches="tight", dpi=300)
plt.close(fig)


# In[ ]:




