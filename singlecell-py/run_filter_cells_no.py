#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

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
# ========== Global Variables ==========

PCs = 20
res = 0.3
top_genes = 2000   # number of top variable genes
dataset = 'S-AG'  # sample name, will affect input/output file names
h5adpath = f'{dataset}_nohomo.h5ad'  # path to read the h5ad file
vars_use = ['sample']  # variables used in Harmony

# List of clusters to keep (example). Adjust these to the cluster IDs you want to retain.
clusterlist = [2,4,5,11]

# Whether to update new cell type annotation
update_cell_type = False
NewCellType = {
    "": [2, 4, 7, 6, 10]
}


# Genes for feature-plot (as in your original script)
Featuregenes = [
    'ESR1', 'EPCAM', 'TOP2A', 'ACTA2', 'ELF5', 'TCF4', 'KRT1', 'PRLR', 'AR', 'PIGR',
    'CD69', 'ADIPOQ', 'LUM', 'PGR', 'VIM', 'PTPRC', 'LALBA', 'LEF1', 'TPM2', 'KRT23',
    'KRT10', 'FAAH', 'TM4SF1', 'PPL', 'WNT11', 'KRTDAP', 'SBSN', 'DSP', 'RAB25',
    'AQP3', 'SHH', 'ATP1A1', 'ATP1B1', 'PROCR', 'WNT6', 'KRT5', 'TP63', 'LGR5', 'KRT5',
    'PIP', 'SOX9', 'SOX10', 'AREG'
]

# Whether to plot feature genes
doFeatureplot = True

# Whether to run Harmony
runharmony = True

random.seed(2024)
np.random.seed(2024)

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

    adata.obs['newcelltype'] = adata.obs[clustering_column].map(cluster_to_celltype)
    print("Finished updating cell annotation:")
    print(adata.obs['newcelltype'].value_counts())
    return adata


def run_harmony(adata, vars_use):
    """
    Perform Harmony batch correction on PCA embeddings.
    """
    print("Running Harmony...")
    pca_result = adata.obsm['X_pca']
    ho = hm.run_harmony(pca_result, adata.obs, vars_use, random_state=42)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print("Harmony finished.")
    return adata


def run_preprocess(adata, n_top_genes=3000):
    """
    Standard Scanpy preprocessing: normalize_total -> log1p -> highly_variable_genes.
    """
    #sc.pp.normalize_total(adata)
    #sc.pp.log1p(adata)
    adata.X=adata.layers["normalized"]
    adata.raw = adata
    #adata.layers["normalized"] = adata.X.copy()
    #print("Data normalization and log1p completed.")
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
    sc.pl.highly_variable_genes(adata)
    print("Highly variable gene selection completed.")


def run_reduce_dimension(adata, run_harmony_flag=False, PCs=20, resolution=0.3, vars_use=None):
    """
    Scale data, PCA, neighbors, Leiden clustering, UMAP.
    Optionally run Harmony if run_harmony_flag is True.
    """
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, use_highly_variable=True)
    sc.pl.pca_variance_ratio(adata, log=False)
    print("PCA completed.")

    if run_harmony_flag and vars_use is not None:
        adata = run_harmony(adata, vars_use=vars_use)
        sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_pcs=PCs, random_state=42, n_neighbors=15)
        print("Neighbors graph constructed using Harmony output.")
    else:
        sc.pp.neighbors(adata, n_pcs=PCs, random_state=42)
        print("Neighbors graph constructed without Harmony.")

    sc.tl.leiden(adata, resolution=resolution, random_state=42)
    print("Leiden clustering done.")
    sc.tl.umap(adata, random_state=42)
    print("UMAP embedding finished.")

    return adata


def do_umap_plots(adata, sample_prefix="sample", feature_genes=None, do_feature=False):
    """
    Plot UMAP with basic annotations and optional feature genes.
    """
    sc.pl.umap(
        adata,
        color=[ "sample", "leiden", "cluster"],
        wspace=0.5,
        size=3,
        ncols=3,
        save=f'{sample_prefix}_nohomo_basic_umap.png'
    )
    if do_feature and feature_genes is not None:
        sc.pl.umap(
            adata,
            color=[*feature_genes,'cluster','sample','leiden'],
            wspace=0.5,
            size=3,
            ncols=4,
            save=f'{sample_prefix}_nohomo_feature_umap.png'
        )
        print("Feature gene UMAP plots completed.")
    else:
        print("Skipping feature gene UMAP plots.")


def export_deg_result(adata):
    """
    Extract DE results from adata.uns['rank_genes_groups'] into a DataFrame.
    Returns only pval < 0.01 subset.
    """
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    data = []
    for group in groups:
        genes = result['names'][group]
        logfoldchanges = result['logfoldchanges'][group]
        pvals = result['pvals'][group]
        pvals_adj = result['pvals_adj'][group]
        for gene, lfc, pval, pval_adj in zip(genes, logfoldchanges, pvals, pvals_adj):
            data.append([group, gene, lfc, pval, pval_adj])
    df = pd.DataFrame(data, columns=['group', 'gene', 'logfoldchange', 'pval', 'pval_adj'])
    df1 = df[df['pval'] < 0.01]
    return df1


def do_deg(adata, sample_prefix="sample"):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        n_genes=5,
        save=f'{sample_prefix}_dotplot_leiden.png',
        min_logfoldchange=0.25
    )
    df1 = export_deg_result(adata)
    df1.to_csv(f'{sample_prefix}_nohomo_ranked_genes_leiden.csv', index=False)
    print("DE analysis completed and results saved.")


def do_cell_barplot(adata, sample_prefix="sample"):
    """
    Plot stacked bar chart of celltype composition by stage.
    """
    if 'newcelltype' not in adata.obs:
        print("Cannot plot stacked bar chart: 'newcelltype' not found in adata.obs.")
        return

    celltype_counts = adata.obs.groupby(['stage', 'newcelltype']).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100

    stage_names = celltype_percentages.index
    cell_types = celltype_percentages.columns

    # Attempt to retrieve color palette from adata; otherwise, fall back to a default palette
    default_palette = adata.uns.get('newcelltype_colors', sns.color_palette("magma", len(cell_types)))

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


# ========== Main Execution ==========

if __name__ == "__main__":
    print(f"Reading h5ad file: {h5adpath}")
    adata_raw = sc.read_h5ad(h5adpath)
    # Keep only specified clusters
    if 'leiden' not in adata_raw.obs:
        print("Warning: 'leiden' was not found in obs. Please check your cluster column.")
    else:
        print(f"Filtering clusters, keeping only: {clusterlist}")
        clusterlist_to_keep=[str(i) for i in clusterlist]
        adata_raw = adata_raw[adata_raw.obs['leiden'].isin(clusterlist_to_keep)].copy()
    adata_raw.obs['cluster'] =adata_raw.obs['leiden']
    # Preprocess
    run_preprocess(adata_raw, n_top_genes=top_genes)

    # Run dimension reduction
    adata_raw = run_reduce_dimension(
        adata_raw,
        run_harmony_flag=runharmony,
        PCs=PCs,
        resolution=res,
        vars_use=vars_use
    )

    # Update celltype annotation if needed
    adata_raw = update_cell_annotations(
        adata_raw,
        clustering_column='leiden',
        new_cell_type_dict=NewCellType,
        update_cell_type=update_cell_type
    )
    output_filename = f'{dataset}_no_test.h5ad'
    adata_raw.write(output_filename)
    # UMAP plots
    do_umap_plots(
        adata=adata_raw,
        sample_prefix=dataset,
        feature_genes=Featuregenes,
        do_feature=doFeatureplot
    )
    sc.tl.dendrogram(adata_raw, groupby='leiden')
    # DE analysis
    do_deg(adata_raw, sample_prefix=dataset)

    # Barplot for new celltype composition
    do_cell_barplot(adata_raw, sample_prefix=dataset)

    # Write out final result
    output_filename = f'{dataset}_no_cleaned.h5ad'
    adata_raw.write(output_filename)
    print(f"All steps completed. The new h5ad file is saved as: {output_filename}")

