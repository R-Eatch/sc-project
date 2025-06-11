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
import louvain
import loompy
print('import successful')


# In[ ]:


adata1 = sc.read_h5ad('../R-MG/1.subset/R-MG_cleaned.h5ad')


# In[ ]:


adata2 = sc.read_h5ad('../R-AG/1.subset/R-AG_cleaned.h5ad')


# In[ ]:


adata_filtered1 = adata1[adata1.obs["time"] <= 4].copy()


# In[ ]:


adata_filtered2 = adata2[adata2.obs["time"] <= 4].copy()


# In[ ]:


adata_filtered1.write('R-MG-E.h5ad')
adata_filtered2.write('R-AG-E.h5ad')

