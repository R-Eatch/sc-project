#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
# Disable hash randomization for reproducibility
os.environ['PYTHONHASHSEED'] = '0'

import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import harmonypy as hm

print("Modules imported successfully.")
sc.settings.figdir = ''


# In[ ]:


PCs = 20 
res = 1.8
dataset = 'S-MG'
h5adpath = f'../{dataset}_cleaned.h5ad'
celltypelist = ["MaSC",
             "MaSC-t2-sg",
             "LumSEC-Lac",
             "LumSEC-Lip",
             "Epi-Krt7",
             'LumHR'
]
random.seed(12345)
np.random.seed(12345)
NewCellType = {
    "StemCells": [0,5,6],
    "Lum0": [2],"LumHR": [1],"LumSEC": [3],"Lum2": [7,8,9],"LumP1": [4]
}
clusterlist = ['4']
Featuregenes = ['Esr1','Epcam','Top2a','Acta2','Zeb1','Tcf4']
update_cell_type = False
doFeatureplot = True
runharmony = False
subset_celltype = True
do_subset_cluster = False
random_subset_cells = False


# In[ ]:


# ========== Functions ==========

def update_cell_annotations(adata, clustering_column, new_cell_type_dict, update_cell_type):
    """
    Update cell annotations based on a dictionary mapping cluster -> new cell type.
    Only performs updates if update_cell_type=True.
    """
    if not update_cell_type:
        print("Cell annotation update is disabled.")
        return adata

    cluster_to_celltype = {}
    for cell_type, clusters in new_cell_type_dict.items():
        for cluster in clusters:
            # if 'leiden' is stored as string in adata.obs, be sure to convert to str
            cluster_to_celltype[str(cluster)] = cell_type

    adata.obs['sub-newcelltype'] = adata.obs[clustering_column].map(cluster_to_celltype)
    print("Finished updating cell annotation:")
    print(adata.obs['sub-newcelltype'].value_counts())
    return adata

def do_umap_plots(adata, sample_prefix="sample", feature_genes=None, do_feature=False):
    """
    Plot UMAP with basic annotations and optional feature genes.
    """
    sc.pl.umap(
        adata,
        color=["newcelltype", "sub-newcelltype", "stage", "sample", "leiden", "gland"],
        wspace=0.5,
        size=3,
        ncols=3,
        save=f'{sample_prefix}_basic_umap.png'
    )
    sc.pl.umap(
        adata,
        color=["sub-newcelltype","leiden"],legend_loc='on data'
        save=f'{sample_prefix}_leiden_umap.png'
    )
    if do_feature and feature_genes is not None:
        sc.pl.umap(
            adata,
            color=[*feature_genes,'newcelltype','sub-newcelltype',"stage",'sample','leiden','gland'],
            wspace=0.5,
            size=3,
            ncols=4,
            save=f'{sample_prefix}_feature_umap.png'
        )
        print("Feature gene UMAP plots completed.")
    else:
        print("Skipping feature gene UMAP plots.")

def do_deg(adata, sample_prefix="sample",groupkey='leiden'):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby=groupkey, method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=groupkey,
        n_genes=5,
        save=f'{sample_prefix}_dotplot_leiden.png',
        min_logfoldchange=0.25
    )
    
    df1 = sc.get.rank_genes_groups_df(adata,group=None,log2fc_min=0.1,pval_cutoff=0.01)
    df1.to_csv(f'{sample_prefix}_ranked_genes_{groupkey}.csv', index=False)
    print("DE analysis completed and results saved.")

def do_cell_barplot(adata, sample_prefix="sample",obsname='newcelltype'):
    """
    Plot stacked bar chart of celltype composition by stage.
    """
    if obsname not in adata.obs:
        print("Cannot plot stacked bar chart: 'newcelltype' not found in adata.obs.")
        return

    celltype_counts = adata.obs.groupby(['stage', obsname]).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100

    stage_names = celltype_percentages.index
    cell_types = celltype_percentages.columns

    # Attempt to retrieve color palette from adata; otherwise, fall back to a default palette
    default_palette = adata.uns.get(f'{obsname}_colors', sns.color_palette("magma", len(cell_types)))

    bottom = pd.Series([0] * len(stage_names), index=stage_names)
    fig, ax = plt.subplots(figsize=(12, 8))

    for cell_type, color in zip(cell_types, default_palette):
        percentages = celltype_percentages[cell_type]
        ax.bar(stage_names, percentages, bottom=bottom, color=color, label=cell_type)
        bottom += percentages

    ax.set_title('Celltype Composition by Stage', fontsize=16)
    ax.set_xlabel('Stage', fontsize=14)
    ax.set_ylabel('Percentage (%)', fontsize=14)
    ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xticks(rotation=45)
    plt.tight_layout()
    fig.savefig(f'{sample_prefix}_cell_bar_plot.png')
    plt.show()
    print("Stacked bar chart for celltype composition saved.")


# In[ ]:


print(f"Reading h5ad file: {h5adpath}")
adata_raw = sc.read_h5ad(h5adpath)
if subset_celltype and 'newcelltype' in adata_raw.obs:
    print("Filtering by specified newcelltype...")
    adata_raw = adata_raw[adata_raw.obs['newcelltype'].isin(celltypelist)].copy()
sc.tl.leiden(adata_raw,resolution = res,random_state=42)
adata_raw.obs['sub-newcelltype']=adata_raw.obs['newcelltype']
# Update celltype annotation if needed
adata_raw = update_cell_annotations(
    adata_raw,
    clustering_column='leiden',
    new_cell_type_dict=NewCellType,
    update_cell_type=update_cell_type
)

# UMAP plots
do_umap_plots(
    adata=adata_raw,
    sample_prefix=dataset,
    feature_genes=Featuregenes,
    do_feature=doFeatureplot
)
sc.tl.dendrogram(adata_raw, groupby='leiden')
# DE analysis
do_deg(adata_raw, sample_prefix=dataset,groupkey='leiden')
do_deg(adata_raw, sample_prefix=dataset,groupkey='sub-newcelltype')
# Barplot for new celltype composition
do_cell_barplot(adata_raw, sample_prefix=dataset,obsname='sub-newcelltype')

# Write out final result
output_filename = f'{dataset}-sub-newcelltype_cleaned.h5ad'
adata_raw.write(output_filename)
print(f"All steps completed. The new h5ad file is saved as: {output_filename}")

