#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as  sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import harmonypy as hm
print('import successful')
sc.settings.n_jobs = 1


# In[4]:

dataset='ALL-AG'
doFeatureplot = True
sc.settings.figdir = ''
Featuregenes = ['Krt5','Acta2','Krt17','Epcam','Pip','Krt7','Krt23','Mif','Melk']

# In[ ]:

adata= sc.read_h5ad('./ALL-AG_cleaned-5_23_final.h5ad')

if doFeatureplot:
    sc.pl.umap(
        adata,
        color=[*Featuregenes,'newcelltype','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_feature.pdf',
    color_map='seismic',
    )
    print('Featureplot Finished')
else:
    print('666666')
if doFeatureplot:
    sc.pl.umap(
        adata,
        color=[*Featuregenes,'newcelltype','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_feature.png',
    color_map='seismic',
    )
    print('Featureplot Finished')
else:
    print('666666')
