#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import anndata
import dynamo as dyn
import scanpy as  sc
import plotly.express as px
import matplotlib.pyplot as plt
import sys
import plotly
dyn.get_all_dependencies_version()
dyn.configuration.set_figure_params(background='white')
print('import successful')


# In[2]:


import scvelo as scv
scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
#scv.set_figure_params('scvelo')  # for beautified visualization


# In[3]:


def run_plotly(adata):
    umap_data = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
    umap_data['cell_ids'] = adata.obs.index 
    umap_data['time'] = adata.obs['time'] 
    fig = px.scatter(umap_data, x='UMAP1', y='UMAP2', hover_data=['cell_ids', 'stage'],width=800,height=450)
    fig.show()
    fig.write_html('interactive_umap_plot.html')


# In[4]:


def run_raw_data(adata):
    adata.X = adata.layers['counts'] ### Reset matrix X to the original matrix
    return(adata)


# In[5]:


def dropcelltype(adata,celltypelist):
    adata=adata[~adata.obs['newcelltype'].isin(celltypelist),:]
    return adata


# In[6]:


###global variable###
dataset = 'S-MAG'
ad_path = f'../1.subset/{dataset}_cleaned.h5ad'
#ad_path2 = f'../2.scvelo/{dataset}_velo_for_dynamo.h5ad'
subset_celltype = False
cores = 8
PCs=10
max_genes=100
use_fixed_points = False
run_preprocess=True
run_animation = False
run_scvelo = False
run_LAPS = False
use_rawdata=True
run_umap=True
use_negbin=False
stem = 'ASC'
cell1 = 'LumSEC-Lip'
cell2 = 'LumSEC-MG-like'
good_fixed_points = [102,227,94]
celltypelist=['Lum-Fibro-like','Lum-Ker-like','LumSEC type3']


# In[7]:


#ad_path="D:/111/R-AG_cleaned.h5ad"
#ad_path2="D:/111/M-MG_velo_for_dynamo.h5ad"


# In[ ]:


dyn.configuration.set_figure_params('', background='white')


# In[8]:


adata = anndata.read(ad_path)


# In[9]:


adata


# In[10]:


if use_rawdata:
    adata=run_raw_data(adata)


# In[11]:


if subset_celltype:
    adata=dropcelltype(adata=adata,celltypelist=celltypelist)


# In[13]:


adata


# In[14]:


if run_preprocess:
    from dynamo.preprocessing import Preprocessor
    preprocessor = Preprocessor()
    preprocessor.preprocess_adata(adata, recipe="monocle")
    #dyn.pp.recipe_monocle(adata)


# In[15]:


if use_negbin:
    dyn.tl.dynamics(adata, model='stochastic', est_method='negbin',cores=cores)
else:
    dyn.tl.dynamics(adata, model='stochastic', cores=cores)

if run_umap:
    dyn.tl.reduceDimension(adata, basis='umap', enforce=False)
dyn.tl.cell_velocities(adata)
dyn.tl.cell_velocities(adata,basis='pca',method='cosine')


# In[16]:


dyn.tl.cell_wise_confidence(adata)


# In[17]:


dyn.vf.VectorField(adata, basis='umap',cores=cores)
dyn.vf.VectorField(adata, basis='pca',cores=cores)


# In[17]:


#dyn.vf.VectorField(adata, basis='umap',cores=cores,dims=PCs)
#dyn.vf.VectorField(adata, basis='pca',cores=cores,dims=PCs)


# In[18]:


dyn.vf.topography(adata, basis='umap')
dyn.vf.topography(adata, basis='pca')


# In[ ]:


scv.pl.velocity_embedding_stream(adata, basis='umap',color='newcelltype',title=f'{dataset}-UMAP-velo',save=f'{dataset}_velo_umap.png')


# In[ ]:


scv.pl.velocity_embedding_stream(adata, basis='pca',color='newcelltype',title=f'{dataset}-PCA-velo',save=f'{dataset}_velo_pca.png')


# In[ ]:


dyn.pl.topography(
    adata,
    color="newcelltype",
    markersize=200,
    basis="umap",
    fps_basis="umap",
    streamline_alpha=1.0,
    save_show_or_return='return',
    figsize=(6,4),
    background='white',color_key=adata.uns['newcelltype_colors']
)
dyn.pl.save_fig(path='./',ext='png',prefix=f'{dataset}-dynamo-umap')


# In[ ]:


dyn.pl.topography(adata,
                  basis='pca',
                  save_show_or_return='return',
                  figsize=(6,4),
                  background='white'
                 )
dyn.pl.save_fig(path='./',ext='png',prefix=f'{dataset}-dynamo-pca')


# In[ ]:


if use_fixed_points:
    
    Xss, ftype = adata.uns['VecFld_umap']['Xss'], adata.uns['VecFld_umap']['ftype']
    adata.uns['VecFld_umap']['Xss'] = Xss[good_fixed_points]
    adata.uns['VecFld_umap']['ftype'] = ftype[good_fixed_points]
else:
    print('do not pick fixed points',flush=True)


# In[ ]:


dyn.pl.topography(
    adata,
    color="newcelltype",
    markersize=200,
    basis="umap",
    fps_basis="umap",
    streamline_alpha=1.0,
    save_show_or_return='return',
    figsize=(6,4),
    background='white',color_key=adata.uns['newcelltype_colors']
)
dyn.pl.save_fig(path='./',ext='png',prefix=f'{dataset}-dynamo-goodpoint-umap')


# In[ ]:


adata


# In[ ]:


#run_plotly(adata)


# In[ ]:


#dyn.pd.fate(adata, "ATCGATGGTGTA",basis='umap')


# In[ ]:


dyn.ext.ddhodge(adata,basis='umap')
dyn.ext.ddhodge(adata,basis='pca')


# In[ ]:


dyn.vf.Potential(adata)


# In[ ]:


adata


# In[ ]:


#dyn.pl.line_integral_conv(adata,basis='umap',save_show_or_return='show')
#dyn.pl.save_fig(path='./',ext='png',prefix=f'{dataset}-integral_conv')


# In[ ]:


if run_animation:
    from matplotlib import animation
    fig, ax = plt.subplots()
    ax = dyn.pl.topography(adata, ax=ax)
    instance = dyn.mv.StreamFuncAnim(adata=adata, ax=ax)
    anim = animation.FuncAnimation(instance.fig, instance.update, init_func=instance.init_background)
    anim.save('animation.mp4', writer='ffmpeg', dpi=300)


# In[ ]:


##scvelo###


# In[ ]:


if run_scvelo:
    
    scv.tl.velocity_graph(adata, vkey='velocity_S', xkey='M_s', n_jobs=4)
    
    scv.tl.velocity_pseudotime(adata, vkey='velocity_S', 
                               #root_key='end',
                               #end_key='root',
                               n_dcs=20,
                               save=f'{dataset}pseudotime.png'
                              )
    
    scv.pl.velocity_embedding_stream(adata, 
                                     color='velocity_S_pseudotime', 
                                     basis='umap', 
                                     cmap='gnuplot',
                                     save=f'{dataset}-scvelo2.png'
                                    )
    
    scv.tl.score_genes_cell_cycle(adata)
    
    scv.pl.scatter(adata, 
                   color_gradients=['S_score', 'G2M_score'], 
                   palette=['green', 'orange'], 
                   smooth=True, perc=[5, 90],
                  save=f'{dataset}score.png')
    
    toptransition_genes = adata.var['use_for_transition'].sort_values(ascending=False).index[:100]
    
    scv.pl.heatmap(adata, 
                   var_names=toptransition_genes, 
                   sortby='velocity_S_pseudotime', 
                   layer='M_s',
                   color_map='viridis',
                   col_color='newcelltype', 
                   palette='viridis', 
                   n_convolve=100,
                   colorbar=True,
                   col_cluster=False, 
                   row_cluster=False,
                   save=f'{dataset}heatmap.png'
                   )
    


# In[ ]:


###LAPS####least action path analyses##############


# In[ ]:


if run_LAPS:
    print('begining running least action path analyses',flush=True)
else:
    sys.exit()


# In[ ]:


adata


# In[ ]:


dyn.pl.streamline_plot(adata, basis="umap", color="newcelltype")
c0 = dyn.tl.select_cell(adata, "newcelltype", stem)
c1 = dyn.tl.select_cell(adata, "newcelltype", cell1)
c2 = dyn.tl.select_cell(adata, "newcelltype", cell2)
c3 = dyn.tl.select_cell(adata, "newcelltype", cell3)
c4 = dyn.tl.select_cell(adata, "newcelltype", cell4)
c5 = dyn.tl.select_cell(adata, "newcelltype", cell5)


# In[ ]:


adata.uns['VecFld_umap']['Xss']


# In[ ]:


from dynamo.tools.utils import nearest_neighbors
fixed_points = adata.uns['VecFld_umap']['Xss']
c0_indices = nearest_neighbors(fixed_points[0], adata.obsm["X_umap"])
c1_indices = nearest_neighbors(fixed_points[1], adata.obsm["X_umap"])
c2_indices = nearest_neighbors(fixed_points[2], adata.obsm["X_umap"])
c3_indices = nearest_neighbors(fixed_points[3], adata.obsm["X_umap"])
c4_indices = nearest_neighbors(fixed_points[4], adata.obsm["X_umap"])
c5_indices = nearest_neighbors(fixed_points[5], adata.obsm["X_umap"])


# In[ ]:


import matplotlib.pyplot as plt
plt.scatter(*adata.obsm["X_umap"].T)
for indices in [
    c0_indices,
    c1_indices,
    c2_indices,
    c3_indices,
    c4_indices,
    c5_indices,
]:
    plt.scatter(*adata[indices[0]].obsm["X_umap"].T)


# In[ ]:


dyn.tl.neighbors(adata, basis="umap", result_prefix="umap")


# In[ ]:


dyn.dynamo_logger.main_silence()
transition_graph = {}
cell_type = [stem, cell1, cell2,cell3,cell4,cell5]
start_cell_indices = [
    c0_indices,
    c1_indices,
    c2_indices,
    c3_indices,
    c4_indices,
    c5_indices,
]
end_cell_indices = start_cell_indices
for i, start in enumerate(start_cell_indices):
    for j, end in enumerate(end_cell_indices):
        if start is not end:
            print(f'        now are : {start}---{end}              ')
            min_lap_t = True if i == 0 else False
            dyn.pd.least_action(
                adata,
                [adata.obs_names[start[0]][0]],
                [adata.obs_names[end[0]][0]],
                basis="umap",
                adj_key="umap_distances",
                min_lap_t= min_lap_t,
                EM_steps=2,
            )
            dyn.pl.least_action(adata, basis="umap")
            plt.savefig('./' + dataset + '_lap_' + cell_type[i]+ "-" + cell_type[j] + '.png')
            lap = dyn.pd.least_action(
                adata,
                [adata.obs_names[start[0]][0]],
                [adata.obs_names[end[0]][0]],
                basis="pca",
                adj_key="cosine_transition_matrix",
                min_lap_t=min_lap_t,
                EM_steps=2,
            )
            genes = adata.var_names[adata.var.use_for_transition]
            max_genes = max_genes
            genes = genes[:max_genes]
            df=pd.DataFrame(genes,columns=['Genes'])
            df.to_csv('./' + dataset + '__transition_genes_' + cell_type[i]+ "-" + cell_type[j] + '.csv', index=False)
            dyn.pl.kinetic_heatmap(
                adata,
                basis="pca",
                mode="lap",
                genes=genes,
                show_colorbar=True,
                project_back_to_high_dim=True,
                figsize=(30, 50),
                
            )
            plt.savefig('./' + dataset + '_heatmap_' + cell_type[i]+ "-" + cell_type[j] + '.png')
            
            # The `GeneTrajectory` class can be used to output trajectories for any set of genes of interest
            gtraj = dyn.pd.GeneTrajectory(adata)
            gtraj.from_pca(lap.X, t=lap.t)
            gtraj.calc_msd()
            ranking = dyn.vf.rank_genes(adata, "traj_msd")

            print(start, "->", end)
            genes = ranking[:5]["all"].to_list()
            arr = gtraj.select_gene(genes)
            
            dyn.pl.multiplot(lambda k: [plt.plot(arr[k, :]), plt.title(genes[k])], np.arange(len(genes)))
            plt.savefig('./' + dataset + '_' + cell_type[i]+ "-" + cell_type[j] + '.png')

            transition_graph[cell_type[i] + "->" + cell_type[j]] = {
                "lap": lap,
                "LAP_umap": adata.uns["LAP_umap"],
                "LAP_pca": adata.uns["LAP_pca"],
                "ranking": ranking,
                "gtraj": gtraj,
            }


# In[ ]:


develope_keys = [f"{stem}->{cell1}",f'{stem}->{cell2}'
                f'{stem}->{cell3}'
f'{stem}->{cell4}'
f'{stem}->{cell5}']
reprogram_keys = [f"{cell1}->{stem}",f'{cell2}->{stem}'
                 f'{cell3}->{stem}'
f'{cell4}->{stem}'
f'{cell5}->{stem}']
all_keys = develope_keys + reprogram_keys


# In[ ]:


from dynamo.plot.utils import map2color
def plotsave_lap(paths):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax = dyn.pl.topography(
        adata, markersize=300, 
        basis="umap", 
        save_show_or_return="return", 
        ax=ax, 
        fps_basis="umap", 
        color="newcelltype", 
        streamline_alpha=1, 
        frontier=False,color_key=adata.uns['newcelltype_colors']
    )
    #ax = ax[0]
    x, y = 0, 1
    # plot paths
    for path in paths:
        lap_dict = transition_graph[path]["LAP_umap"]
        for prediction, action in zip(lap_dict["prediction"], lap_dict["action"]):
            ax.scatter(*prediction[:, [x, y]].T, c=map2color(action))
            ax.plot(*prediction[:, [x, y]].T, c="k")
    plt.savefig(f'{dataset}-all_LAP.png')


# In[ ]:


plotsave_lap(all_keys)


# In[ ]:


from dynamo.plot.utils import map2color
def plot_lap(paths):
    fig, ax = plt.subplots(figsize=(5, 4))
    ax = dyn.pl.streamline_plot(
        adata, basis="umap", save_show_or_return="return", ax=ax, color="newcelltype", frontier=True,color_key=adata.uns['newcelltype_colors']
    )
    ax = ax[0]
    x, y = 0, 1

    # plot paths
    for path in paths:
        lap_dict = transition_graph[path]["LAP_umap"]
        for prediction, action in zip(lap_dict["prediction"], lap_dict["action"]):
            ax.scatter(*prediction[:, [x, y]].T, c=map2color(action))
            ax.plot(*prediction[:, [x, y]].T, c="k")
    plt.savefig(f'{dataset}-develop_LAP.png')


# In[ ]:


plot_lap(develope_keys)


# In[ ]:


###test code do not run###


# In[ ]:


#adata.obs['leiden']


# In[ ]:


#adata = adata[adata.obs['leiden'].isin(['0','1','2']),:].copy()


# In[ ]:


#adata.var


# In[ ]:


#type(adata.uns['VecFld_umap']['Xss'])


# In[ ]:


#adata.uns['VecFld_umap']['Xss'].type()


# In[ ]:




