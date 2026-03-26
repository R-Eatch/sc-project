#!/usr/bin/env python
# coding: utf-8

# In[1]:


import loompy
import scvelo as scv
import scanpy as sc
import os
import re
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


# In[2]:


scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[ ]:


###global variable###

#mapping_dict = {'S-MG-3MTH': 'S-3MTh-MG-v','S-MG-8M-3': 'S-8M-MG-3-v','S-MG-8M-4': 'S-8M-MG-4-v','S-MG-GES14': 'S-GES14-MG','S-MG-LA14': 'S-LA14-MG','S-MG-P20': 'S-P20-MG-v','S-MG-P7': 'S-P7-MG-v'}
#mapping_dict = { 'S-AG-1MTH-1': 'SG-AG-1Mth-1-v','S-AG-1MTH-2': 'SG-AG-1Mth-2-v','S-AG-8M-2': 'S-8m-AG-2-v','S-AG-8M-3': 'S-8m-AG-3-v','S-AG-GES14': 'S-GES14-AG', 'S-AG-LA-1': 'SG-LA-AG-1-v','S-AG-LA-2': 'SG-LA-AG-2-v', 'S-AG-P20': 'S-P20-AG-v'}
#mapping_dict = { 'M-MG-3WK-1': 'MG-3WK-1', 'M-MG-3WK-2': 'MG-3WK-2','M-MG-8WK-1': 'MG-8WK-1','M-MG-8WK-2': 'MG-8WK-2','M-MG-E13_5': 'M-E13_5-MG','M-MG-E16_5': 'M-E16-5-MG','M-MG-GES13_5': 'M-GES13_5-MG','M-MG-GES16_5': 'M-GES16-5-SG','M-MG-LA-1': 'M-LA-MG-1','M-MG-LA-2': 'M-LA-MG-2', 'M-MG-P1': 'M-P1-MG'}
#mapping_dict ={'R-AG-10WK-1': 'R-10wk-AG-1', 'R-AG-25WK-1': 'R-25wk-AG-1', 'R-AG-25WK-2': 'R-25WK-AG-2-r-bam', 'R-AG-8WK-1': 'R-8wk-AG-1', 'R-AG-E26': 'R-E26-AG', 'R-AG-GES12': 'R-GES12-AG', 'R-AG-GES17': 'R-GES17-AG', 'R-AG-GES23': 'R-GES23-AG', 'R-AG-LA-2': 'R-LA-AG-2', 'R-AG-LA': 'J103660_R-La-AG', 'R-AG-P1': 'R-P1-AG'}
mapping_dict ={'R-MG-23WK-3': 'J103660_R-MG-23WK-3', 'R-MG-23WK-4': 'J103660_R-MG-23WK-4', 'R-MG-8WK-1': 'J103660_R-MG-8WK-1', 'R-MG-8WK-2': 'J103660_R-MG-8WK-2', 'R-MG-E17': 'R-E17-MG', 'R-MG-E23': 'R-E23-MG', 'R-MG-GES12': 'R-GES12-MG', 'R-MG-GES17': 'R-GES17-MG', 'R-MG-GES23': 'R-GES23-MG-r-bam', 'R-MG-LA-2': 'R-LA-MG-2', 'R-MG-LA': 'J103660_R-La-MG', 'R-MG-P1': 'R-P1-MG'}
#mapping_dict ={'R-CG-13WK': 'R-13WK-CG', 'R-CG-23WK': 'R-23WK-CG', 'R-CG-E26': 'R-E26-CG', 'R-CG-GES17': 'R-GES17-CG', 'R-CG-GES30': 'R-GES30-CG', 'R-CG-LA': 'R-LAC-CG', 'R-CG-P1': 'R-P1-CG'}
dataset = 'R-MG'
fill_na = True ##S-MG,S-AG =TRUE
run_geneconvert=True## R,S =TRUE

output_loom = f'{dataset}_loom_merged.loom'
adata_path = f'../{dataset}.h5ad'


# In[ ]:


def gene_convert(ldata,df):
    df1=df[(['sg','Gene name'])]
    df1=df1.rename(columns={'sg':'Accession'})
    df2 = pd.DataFrame(ldata.var['Accession'])
    df2=df2.reset_index()
    df3=pd.merge(df1,df2,on='Accession')
    gene_dict=dict(zip(df3['Gene'],df3['Gene name']))
    print(df3.shape)
    var_names = pd.Series(ldata.var_names, index=ldata.var_names)
    mapped_var_names = var_names.map(gene_dict).fillna(var_names)
    ldata.var_names = mapped_var_names.values
    # 提取Accession列
    set1 = set(df1['Accession'])
    set2 = set(df2['Accession'])
    intersection = set1 & set2  # 重合部分
    only_in_df1 = set1 - set2   # 仅在df1中
    only_in_df2 = set2 - set1   # 仅在df2中
    result = {
        "Intersection": len(list(intersection)),
        "Only in DF1": len(list(only_in_df1)),
        "Only in DF2": len(list(only_in_df2))
    }
    
    result_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in result.items()]))
    print(result)
    return ldata


# In[ ]:


def fill_ldata_with_nan(adata, ldata):
    """
    直接重建 ldata，使其基因与 adata 的基因集合完全对齐，缺失部分填充 NaN，并保持 X 和 layers 为稀疏矩阵。
    
    参数：
        adata (AnnData): 包含完整基因集合的 AnnData 对象。
        ldata (AnnData): 需要调整基因集合的 AnnData 对象。
    
    返回：
        AnnData: 重建后的 ldata 对象，基因与 adata.var_names 对齐，缺失值填充为 NaN。
    """
    # 确保 var_names 唯一
    ldata.var_names_make_unique()
    adata.var_names_make_unique()

    # 创建一个完整的基因集合
    all_genes = pd.Index(adata.var_names).union(pd.Index(ldata.var_names)).unique()

    # 创建新的 layers
    new_layers = {}
    for layer in ldata.layers.keys():
        # 初始化新的稀疏矩阵，形状为 (n_obs, len(all_genes))
        new_layer = csr_matrix((ldata.n_obs, len(all_genes)), dtype=np.float32)

        # 找到 ldata 当前基因的位置索引
        existing_genes_idx = [list(all_genes).index(g) for g in ldata.var_names if g in all_genes]

        # 将已有数据填入稀疏矩阵
        new_layer[:, existing_genes_idx] = csr_matrix(ldata.layers[layer])

        # 添加到新的 layers
        new_layers[layer] = new_layer

    # 创建新的 var 数据框
    new_var = pd.DataFrame(index=all_genes)
    for col in ldata.var.columns:
        new_var[col] = ldata.var[col].reindex(all_genes)

    # 创建新的 ldata 对象，初始化为稀疏矩阵
    new_ldata = sc.AnnData(
        X=csr_matrix((ldata.n_obs, len(all_genes))),  # 初始化为空的稀疏矩阵
        obs=ldata.obs.copy(),
        var=new_var
    )

    # 添加稀疏 layers
    for layer, data in new_layers.items():
        new_ldata.layers[layer] = data

    return new_ldata


# In[ ]:


### merge loom file ###
loom_files = []
for key,value in mapping_dict.items():
    loom_filename = f'../../data/loom/{key}.loom'
    #loom_filename = f'D:/BaiduNetdiskDownload/loom/{key}.loom'
    loom_files.append(loom_filename)
if os.path.exists(output_loom) and os.path.getsize(output_loom) > 0:
    print('loom exists and is valid, skip!',flush=True)
else:
    print(f'merging loom file of {dataset} {"#"*50}',flush=True)
    loompy.combine(loom_files, output_loom,key="Accession")
    print(f'merge loom file completed {dataset}',flush=True)


# In[ ]:


ldata = scv.read(output_loom)
adata = sc.read_h5ad(adata_path)
print("here are the stupid loom data")
print(ldata)


# In[ ]:


if run_geneconvert:
    print('convert human genes to mouse genes')
    ldata = gene_convert(ldata=ldata,df=pd.read_csv('../../data/ortho_mrs.csv'))#ortho_mrs.csv


# In[ ]:


#ldata = scv.read('D:/111/M-MG_loom_merged.loom')
#adata = sc.read_h5ad('D:/111/M-MG.h5ad')


# In[ ]:


###  Standardize adata cell ID format ###
adata.obs['cellid'] = (
    adata.obs['sample'].astype(str) + '-' + adata.obs_names.astype(str)
).str.replace(r'(-\d+)+$', '', regex=True)
len(adata.obs['cellid'].unique())


# In[ ]:


### Standardize loomfile cell ID format ###
loom_cell_ids = ldata.obs_names
new_loom_cell_ids = loom_cell_ids.copy()
for key, value in mapping_dict.items():
    new_loom_cell_ids = [cell_id.replace(value, key).replace(":", "-").rstrip('x') for cell_id in new_loom_cell_ids]
ldata.obs_names = new_loom_cell_ids


# In[ ]:


if fill_na:
    ldata=fill_ldata_with_nan(adata=adata,ldata=ldata)


# In[ ]:


print(adata)


# In[ ]:


### merge adata and loom ###
h5ad_uids = set(adata.obs['cellid'])
ldata_cellid = set(ldata.obs_names)
common_cells = h5ad_uids.intersection(ldata_cellid)
if len(common_cells) != len(h5ad_uids):
    print(f"WARNING:The number of common cells ({len(common_cells)}) does not match the number of cells in 'adata' ({len(h5ad_uids)}).")
common_genes = adata.var_names.intersection(ldata.var_names)
adata_sub = adata[adata.obs['cellid'].isin(common_cells),adata.var_names.isin(common_genes)]
print
print(adata_sub)
ldata_sub = ldata[ldata.obs_names.isin(common_cells),ldata.var_names.isin(common_genes)]
print(ldata_sub)
adata_sub.obs_names=adata_sub.obs['cellid']
print(f'{len(common_cells)},{len(common_genes)} {"#"*80}',flush=True)
adata1 = scv.utils.merge(adata_sub, ldata_sub)
print(adata1)
adata1.write(f'{dataset}_for_subset.h5ad')


# In[ ]:


###test code do not run ##########


# In[ ]:


#################################

