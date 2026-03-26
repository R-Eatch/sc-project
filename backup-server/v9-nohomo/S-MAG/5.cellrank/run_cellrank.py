#!/usr/bin/env python
"""
Run CellRank on sugar‑glider MG/AG epithelium
============================================
This script extends the original **run_cellrank.py** by first importing
pseudotime values that were computed in Monocle3 and saved as a CSV
(`S_MAG_pseudotime.csv`), writing them into `adata_pala.obs`, and then
proceeding with the standard CellRank workflow (VelocityKernel +
PseudotimeKernel, GPCCA, fate probabilities, driver genes, etc.).

Usage (bash):
    python run_cellrank_S-MAG.py 

Make sure the following files are present (relative paths can be changed
below):
    •  ../6.monocle3/S_MAG.pseudotime.h5ad      – AnnData with UMAP etc.  
    •  ../2.scvelo/S-MG_velofinished.h5ad       – AnnData with velocity   
    •  S_MAG_pseudotime.csv                     – two columns: cell,pseudotime

Author: Adrian (generated via ChatGPT‑o3)
"""

# -----------------------------------------------------------------------------
# 0 ‑‑ imports & settings
# -----------------------------------------------------------------------------
import pandas as pd
import scanpy as sc
import cellrank as cr
import scvelo as scv
import numpy as np
import warnings

warnings.simplefilter("ignore", category=UserWarning)

# nice plotting defaults
sc.settings.set_figure_params(frameon=False, dpi=200)
scv.settings.presenter_view = True
cr.settings.verbosity = 2

# -----------------------------------------------------------------------------
# 1 ‑‑ global params (edit here)                                                
# -----------------------------------------------------------------------------
dataset          = "S-MAG"          # shorthand used in file names
cores            = 8               # parallel threads
pseudotime_csv   = "S_MAG_pseudotime.csv"  # path to Monocle3 pseudotime CSV

terminal_states  = ['LumSEC-MG-like', 'LumSEC-Lip', 'LumSEC', 'LumHR', 'Epi-SBSN', 'Basal']
initial_state    = "StemCells"

n_genes          = 50              # top lineage‑driver genes per lineage
pseudotime_key   = "pseudotime"    # column name to be added to adata.obs
n_states         = [6, 12]         # GPCCA macrostates (coarse, fine)

cr.settings.n_jobs = cores

# paths to AnnData files -------------------------------------------------------
adata_velo_path  = f"../2.scvelo/{dataset}_velofinished.h5ad"
#adata_pala_path  = f"../6.monocle3/{dataset}.pseudotime.h5ad"

# -----------------------------------------------------------------------------
# 2 ‑‑ load AnnData & inject pseudotime                                         
# -----------------------------------------------------------------------------
print("▶ Loading AnnData objects …")
adata_velo = sc.read_h5ad(adata_velo_path)
adata_pala = adata_velo

print(f"   • adata_velo: {adata_velo.n_obs} cells, {adata_velo.n_vars} genes")
print(f"   • adata_pala: {adata_pala.n_obs} cells, {adata_pala.n_vars} genes\n")

print("▶ Importing Monocle3 pseudotime from CSV …")
pt_df = pd.read_csv(pseudotime_csv)
if not {"cell", "pseudotime"}.issubset(pt_df.columns):
    raise ValueError("CSV must contain 'cell' and 'pseudotime' columns.")
pt_df.set_index("cell", inplace=True)

# align & write
after = adata_pala.obs_names.intersection(pt_df.index)
print(f"   Matched {len(after)} / {adata_pala.n_obs} cells.")
adata_pala.obs[pseudotime_key] = np.nan  # initialise column with NaNs
adata_pala.obs.loc[after, pseudotime_key] = pt_df.loc[after, "pseudotime"].values

if adata_pala.obs[pseudotime_key].isna().all():
    raise RuntimeError("No overlapping barcodes between AnnData and CSV.")

print("   Pseudotime column added to adata_pala.obs as '" + pseudotime_key + "'.\n")

# -----------------------------------------------------------------------------
# 3 ‑‑ build kernels                                                            
# -----------------------------------------------------------------------------
print("▶ Computing VelocityKernel …")
scv.logging.print_version()
scv.settings.verbosity = 3
vk = cr.kernels.VelocityKernel(adata_velo, n_jobs=cores)
vk.compute_transition_matrix(n_jobs=cores)
print(vk, "\n")

print("▶ Computing PseudotimeKernel …")
pk = cr.kernels.PseudotimeKernel(
    adata_pala, time_key=pseudotime_key, n_jobs=cores
)
pk.compute_transition_matrix(n_jobs=cores)
print(pk, "\n")

# -----------------------------------------------------------------------------
# 4 ‑‑ combined kernel & GPCCA                                                  
# -----------------------------------------------------------------------------
cbk = 0.5 * pk + 0.5 * vk

print("▶ Fitting GPCCA estimator …")
g = cr.estimators.GPCCA(cbk)
print(g)

g.fit(cluster_key="newcelltype", n_states=n_states)

# plot macrostates -------------------------------------------------------------
print("▶ Plotting macrostates …")
cr.settings.figdir = ""
sc.settings.figdir = ""

g.plot_macrostates(which="all", legend_loc="right", s=100,
                   save=f"{dataset}-all-macrostates.png")

g.set_terminal_states(states=terminal_states)
g.plot_macrostates(which="terminal", legend_loc="right", s=100,
                   save=f"{dataset}-terminal-macrostates.png")

g.set_initial_states(states=initial_state)
g.plot_macrostates(which="initial", legend_loc="right", s=100,
                   save=f"{dataset}-initial-macrostates.png")

# -----------------------------------------------------------------------------
# 5 ‑‑ fate probabilities & lineage drivers                                     
# -----------------------------------------------------------------------------
print("▶ Computing fate probabilities …")
g.compute_fate_probabilities(n_jobs=cores)
g.plot_fate_probabilities(same_plot=False,
                          save=f"{dataset}-fate_probabilities.png")

g.plot_fate_probabilities(same_plot=True,
                          save=f"{dataset}-celltype-fate_probabilities.png")

# lineage driver genes & GAM heatmaps -----------------------------------------
model = cr.models.GAM(adata_pala, n_knots=cores)

for cellname in terminal_states:
    print(f"▶ Processing lineage: {cellname} …")
    delta_df = g.compute_lineage_drivers(lineages=[cellname],
                                         cluster_key="newcelltype",
                                         n_jobs=cores)
    delta_df.to_csv(f"{cellname}-driver-gene.csv")

    adata_pala.obs[f"fate_prob_{cellname}"] = g.fate_probabilities[cellname].X.flatten()

    sc.pl.embedding(
        adata_pala,
        basis="umap",
        color=[f"fate_prob_{cellname}"] + list(delta_df.index[:8]),
        color_map="viridis", s=50, ncols=3,
        vmax="p96", save=f"{dataset}-{cellname}-driver-genes.png")

    cr.pl.heatmap(
        adata_pala,
        model=model,
        lineages=cellname,
        data_key="MAGIC_imputed_data",  # ensure this layer exists
        genes=delta_df.head(n_genes).index,
        time_key=pseudotime_key,
        figsize=(12, 10),
        show_all_genes=True,
        weight_threshold=(1e-3, 1e-3),
        n_jobs=cores,
        save=f"{dataset}-{cellname}-cr-heatmap.png")

print("✔ All done – results saved to current directory.")

