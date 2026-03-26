#!/usr/bin/env python
# coding: utf-8

# In[2]:


from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix


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


def clean_h5ad(dataset,df,sp,obslist):
    adata=sc.read_h5ad(f"/data01/sunxuebo/project/scrnaseq/v8-python/{dataset}/1.subset/{dataset}_cleaned.h5ad")
    if sp == 'ms':
        mapping=dict(zip(df["Gene name"],df["Gene stable ID"]))
    if sp =='rb':
        mapping=dict(zip(df["Gene name"],df["Rabbit gene stable ID"]))
    if sp=='sg':
        mapping=dict(zip(df["Gene name"],df["sg"]))
    adata.X=adata.layers['normalized']
    adata.X = csr_matrix(adata.X)
    adata.var_names = adata.var_names.map(mapping)
    del adata.layers
    adata=dropobs(adata,obslist,sp=sp)
    adata.write(f"./{dataset}_counts.h5ad")
    print(adata)
    return None   


# In[ ]:


###########global variable############


# In[23]:


obslist=['newcelltype','celltype','gland','species']
cut_df=False ## first run set true
sp1='rb'
sp2='sg'
gd1='MG'
gd2='MG'
dataset1='R-MG'
dataset2='S-MG'


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


# In[64]:


clean_h5ad(dataset=dataset1,df=df,sp=sp1,obslist=obslist)


# In[65]:


clean_h5ad(dataset=dataset2,df=df,sp=sp2,obslist=obslist)


# In[ ]:


#####SAMAP-generate SAM file#############


# In[75]:


resolutions={sp1:0.3,sp2:0.3}


# In[76]:


fn1 = f'./{dataset1}_counts.h5ad'
fn2 = f'./{dataset2}_counts.h5ad'


# In[77]:


filenames = {sp1:fn1,sp2:fn2}


# In[78]:


sams = {sp1:fn1,sp2:fn2}
sm = SAMAP(
        sams,
        f_maps = '/data01/sunxuebo/project/scrnaseq/v8-python/samap/re-maps/' ,save_processed=True,resolutions=resolutions
    )
###please use server to run SAMAP###

