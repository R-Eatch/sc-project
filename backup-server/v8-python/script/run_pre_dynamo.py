#!/usr/bin/env python
# coding: utf-8

# In[1]:


import loompy
import scvelo as scv
import scanpy as sc
import os


# In[2]:


scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[3]:


def run_raw_data(adata):
    import anndata
    raw_data = adata.raw.X 
    raw_data = raw_data if isinstance(raw_data, np.ndarray) else raw_data.toarray()  
    new_adata = anndata.AnnData(
        X=raw_data,  
        obs=adata.obs,  
        var=adata.raw.var  
    )
    if 'X_umap' in adata.obsm:
        new_adata.obsm['umap'] = adata.obsm['X_umap']
    return(new_adata)


# In[ ]:


###global variable###
datasets = ['M-E13-5-MG',
            'M-E16-5-MG',
            'MG-3WK-1',
            'MG-3WK-2',
            'MG-8WK-1',
            'MG-8WK-2',
            'M-GES13-5-MG',
            'M-GES16-5-MG',
            'M-LA-MG-1',
            'M-LA-MG-2',
            'M-P1-MG']
mapping_dict = {
    'M-MG-8WK-2': 'MG-8WK-2',
    'M-MG-3WK-2': 'MG-3WK-2',
    'M-MG-LA-1': 'M-LA-MG-1',
    'M-MG-P1': 'M-P1-MG',
    'M-MG-GES16-5': 'M-GES16-5-SG',###This is not a misprint,this dataset is MG###
    'M-MG-E13_5': 'M-E13_5-MG',
    'M-MG-3WK-1': 'MG-3WK-1',
    'M-MG-8WK-1': 'MG-8WK-1',
    'M-MG-E16-5': 'M-E16-5-MG',
    'M-MG-GES13_5': 'M-GES13_5-MG',
    'M-MG-LA-2': 'M-LA-MG-2'
}
dataset = 'M-MG'
output_loom = f'{dataset}_loom_merged.loom'
loom_path = '../data/'
#adata_path = f'../{dataset}.h5ad'
ad_path = f'../1.subset/{dataset}_cleaned.h5ad'


# In[ ]:


loom_files = []
for name in datasets:
    loom_filename = f'{name}.loom'
    loom_files.append(loom_filename)
loom_files


# In[ ]:


if os.path.exists(output_loom) and os.path.getsize(output_loom) > 0:
    print('loom exists and is valid, skip!',flush=True)
else:
    print(f'merging loom file of {dataset} {"#"*50}',flush=True)
    loompy.combine(loom_files, output_loom)
    print(f'merge loom file completed {dataset}',flush=True)


# In[ ]:


ldata = scv.read(output_loom)
ad = sc.read_h5ad(ad_path)
adata=run_raw_data(ad)


# In[ ]:


loom_cell_ids = ldata.obs_names
new_loom_cell_ids = loom_cell_ids.copy()
for key, value in mapping_dict.items():
    new_loom_cell_ids = [cell_id.replace(value, key).replace(":", "__").rstrip('x') for cell_id in new_loom_cell_ids]
ldata.obs_names = new_loom_cell_ids


# In[ ]:


h5ad_uids = set(adata.obs['uid'])
adata.obs_names = adata.obs['uid']
ldata_cellid = set(ldata.obs_names)
common_cells = h5ad_uids.intersection(ldata_cellid)
common_genes = adata.var_names.intersection(ldata.var_names)
adata_sub = adata[adata.obs['uid'].isin(common_cells),adata.var_names.isin(common_genes)].copy()
ldata_sub = ldata[ldata.obs_names.isin(common_cells),ldata.var_names.isin(common_genes)].copy()
print(f'{len(common_cells)},{len(common_genes)} {"#"*80}',flush=True)
adata1 = scv.utils.merge(adata_sub, ldata_sub)
adata1.write(f'{dataset}_for_dynamo.h5ad')

