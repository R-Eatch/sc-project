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
h5adpath = f'S-MG-sub-newcelltype_cleaned.h5ad'
celltypelist = [
             "Epi-Krt7",
             #'LumHR'
]
random.seed(12345)
np.random.seed(12345)
NewCellType = {
    'Epi-Krt7-c1': [19],
    'Epi-Krt7-c2': [14],
    'Epi-Krt7-c3': [17],
    'Epi-Krt7-c4': [7],
    'Epi-Krt7-c5': [20],

    'LumSEC-Lac-c1': [9],
    'LumSEC-Lac-c2': [3],
    'LumSEC-Lac-c3': [4],
    'LumSEC-Lac-c4': [12],
    'LumSEC-Lac-c5': [11],
    'LumSEC-Lac-c6': [18],
    'LumSEC-Lac-c7': [28],
    'LumSEC-Lac-c8': [24],

    'LumSEC-Lip-c1': [16],
    'LumSEC-Lip-c2': [23],

    'LumHR-c1': [13],
    'LumHR-c2': [5],
    'LumHR-c3': [10],
    'LumHR-c4': [2],
    'LumHR-c5': [8],
    'LumHR-c6': [6],
    'LumHR-c7': [25],
    'LumHR-c8': [22],
    'LumHR-c9': [26],
    'LumHR-c10': [30],

    'MaSC-c1': [0],
    'MaSC-c2': [1],
    'MaSC-c3': [15],
    'MaSC-c4': [27],
    'MaSC-c5': [21],
    'MaSC-c6': [29],


}
clusterlist = ['4']
Featuregenes = ['Esr1','Epcam','Top2a','Acta2','Zeb1','Tcf4']
update_cell_type = True
doFeatureplot = True
runharmony = False
subset_celltype = True
do_subset_cluster = False
random_subset_cells = False


# In[ ]:


# ========== Functions ==========


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
    fig.savefig(f'{sample_prefix}_cell_bar_plot-KRT7.png')
    plt.show()
    print("Stacked bar chart for celltype composition saved.")


# In[ ]:


print(f"Reading h5ad file: {h5adpath}")
adata_raw = sc.read_h5ad(h5adpath)
if subset_celltype and 'newcelltype' in adata_raw.obs:
    print("Filtering by specified newcelltype...")
    adata_raw = adata_raw[adata_raw.obs['newcelltype'].isin(celltypelist)].copy()

# Barplot for new celltype composition
do_cell_barplot(adata_raw, sample_prefix=dataset,obsname='sub-newcelltype')

#print(f"All steps completed. The new h5ad file is saved as: {output_filename}")

