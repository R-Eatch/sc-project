import cell2fate as c2f
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import os


###global variable###
dataset = 'S-MG'
ad_path = f'../1.subset/{dataset}_cleaned.h5ad'

adata = sc.read_h5ad(ad_path)
adata.X=adata.layers['counts']
clusters_to_remove = []
adata =  c2f.utils.get_training_data(adata, cells_per_cluster = 10**5, cluster_column = 'leiden',
                                    remove_clusters = clusters_to_remove,
                                min_shared_counts = 20, n_var_genes= 3000)
c2f.Cell2fate_DynamicalModel.setup_anndata(adata, spliced_label='spliced', unspliced_label='unspliced')
n_modules = c2f.utils.get_max_modules(adata)

mod = c2f.Cell2fate_DynamicalModel(adata, n_modules = n_modules)

mod.train()

adata = mod.export_posterior(adata)

adata = mod.compute_module_summary_statistics(adata)

mod.plot_module_summary_statistics(adata, save = 'module_summary_stats_plot.pdf')

mod.compare_module_activation(adata, chosen_modules = [1,2,3,4],
                         save = 'module_activation_comparison.pdf')

mod.compute_and_plot_total_velocity_scvelo(adata, save = 'total_velocity_plots.png', delete = False)

