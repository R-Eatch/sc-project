#!/usr/bin/env python
# coding: utf-8

# In[22]:


import scanpy as sc
import gc
##############################
genes_to_plot = [
    "Sbsn", "Dmkn", "Rbp2", "Serpinb2", "Sult2b1", "Cyp4f39", "Ovol1",
    "Ly6g6c", "Cers3", "Gltp", "Cryba4", "Spink5", "Scin", "Iffo2",
    "Tmem40", "Il1rl2", "Clip4", "Dsp", "Ggh", "Prom2", "Rhov", "Cltb",
    "Hopx", "Mall", "Mboat1", "Lrrc1", "Trim29", "Pkp1"
]
obs_list=['newcelltype','stage']
##############################
datasetlist = ['M-MG',
               'R-MG',
               'S-MG', 
               'R-AG',
                'R-CG',
                'S-AG'
              ]
data_path = "/data01/sunxuebo/project/scrnaseq/v8-python/"


# In[ ]:


for dataset in datasetlist:
    print(f"正在处理数据集：{dataset}")
    adata = sc.read_h5ad(f"{data_path}{dataset}/1.subset/{dataset}_cleaned.h5ad")
    adata.X=adata.layers['normalized']
    valid_genes = [gene for gene in genes_to_plot if gene in adata.var_names or gene == 'newcelltype' or gene == 'stage']
    valid_genes_for_umap=valid_genes+obs_list
    sc.pl.umap(
        adata,
        color=valid_genes_for_umap,
        color_map='viridis',
        ncols=3,
        frameon=False,
        legend_loc='on data',
        save=f"_{dataset}_featureplot.png",
        show=False
    )
    # sc.pl.embedding(
    # adata,
    # basis='umap',         # Specify which embedding to use (essential!)
    # color='pseudotime',
    # color_map='plasma',
    # frameon=False,
    # colorbar_loc=None,
    # title=f'{dataset}_pseudotime',# This is the key parameter to remove the colorbar
    # save=f"_{dataset}_pseudotime-umap.pdf", # Consider a new filename
    # show=True
    # # You might experiment with the 'size' parameter if points are too large/small
    # # size=10 # Example: Adjust marker size (default depends on data size)
# )
    for obs1 in obs_list:
        sc.pl.violin(
        adata,
        keys= valid_genes,
        groupby=obs1, 
        jitter=0.4,
        rotation=90,
        size=2,stripplot=False,save=f"_{dataset}_{obs1}_violin.png"
        )
    del adata
    gc.collect()
print("全部绘制完成，带标签的图片已保存至figures目录中。")

