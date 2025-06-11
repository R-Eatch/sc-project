#!/usr/bin/env python
# coding: utf-8

# In[4]:


import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import scipy
import scanpy as  sc
import anndata as ad
import warnings
import seaborn as sns
warnings.filterwarnings("ignore", category=UserWarning, message="ImplicitModificationWarning")
print('import successful')


# In[ ]:


####global variable####
path="../../h5ad2seurat/allspecies_v8_nosg_raw_re2.h5ad"


# In[1]:


def draw_umap_plot(adata,dataset):
    sc.settings.figdir = ''
    sc.pl.umap(
        adata,
        color=["celltype", "species", "stage",'sample','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_ac_celltype.png',
    color_map='viridis'
    )
    plt.close()


# In[2]:


def do_DEG(ea,dataset):
    sc.tl.rank_genes_groups(ea,groupby='celltype',method = 'wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ea,groupby='celltype',n_genes=5,save=f'{dataset}_ac_dotplot.png',min_logfoldchange=0.25)
    plt.close()


# In[3]:


def do_cell_barplot(ea,dataset):
    celltype_counts = ea.obs.groupby(['stage', 'celltype']).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    celltype_percentages_long = celltype_percentages.reset_index().melt(id_vars='stage', var_name='Celltype', value_name='Percentage')
    default_palette = ea.uns['newcelltype_colors']
    plt.figure(figsize=(14, 8))
    sns.set(style="whitegrid")
    
    times = celltype_percentages.index
    cell_types = celltype_percentages.columns
    
    colors = sns.color_palette("magma", len(cell_types))
    colors = default_palette
    
    bottom = pd.Series([0] * len(times), index=times)
    
    
    fig3, ax = plt.subplots(figsize=(14, 8))
    
    for cell_type, color in zip(cell_types, colors):
        percentages = celltype_percentages[cell_type]
        ax.bar(times, percentages, bottom=bottom, color=color, label=cell_type)
        bottom += percentages
    
    ax.set_title(f'Percentage of Each Celltype Over Time - {dataset}', fontsize=16)
    ax.set_xlabel('Stage', fontsize=14)
    ax.set_ylabel('Percentage', fontsize=14)
    
    ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    fig3.savefig(f"{dataset}_ac_bar_plot.png")
    plt.close()


# In[ ]:


def split_h5ad_file(input_file):
    adata = sc.read_h5ad(input_file)
    conditions = {
        "M-MG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "MG"),
        "R-MG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "MG"),
        "R-AG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "AG"),
        "S-MG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "MG"),
        "S-AG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "AG"),
        "R-CG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "CG")
    }
    sc.settings.figdir = ''
    for i, (group_name, condition) in enumerate(conditions.items()):
        subset_adata = adata[condition].copy()
        # 在每个子图中绘制 UMAP 图
        sc.pl.umap(subset_adata,
                   color=["celltype"],
                   title=group_name,  
                   show=True
                  save=f'{group_name}_ac1.png')  

        # 调用其他函数并保存图像，但不影响当前子图
        draw_umap_plot(subset_adata, group_name)  # 其他绘图函数会保存自己的图像
        do_DEG(subset_adata, group_name)          # DEG 绘图函数会保存自己的图像
        do_cell_barplot(subset_adata, group_name) # 条形图绘图函数会保存自己的图像

    plt.close()


# In[ ]:


if __name__ == "__main__":
    split_h5ad_file(path)

