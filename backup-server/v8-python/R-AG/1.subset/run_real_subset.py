#!/usr/bin/env python
# coding: utf-8

# In[110]:


import os
# Set the environment variable to disable hash randomization
os.environ['PYTHONHASHSEED'] = '0'
#If you find discrepancies between the UMAP and Leiden results, indicating reproducibility issues, 
#please refer to this link: https://github.com/scverse/scanpy/issues/1009 for potential solutions.
#https://umap-learn.readthedocs.io/en/latest/reproducibility.html
# set export OMP_NUM_THREADS=1 Multithreading can lead to reproducibility issues!!!
#Theoretically, minor differences in UMAP and inconsistencies in Leiden clustering do not affect their biological significance.
import scanpy as  sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import loompy
import harmonypy as hm
print('import successful')


# In[111]:


###global variable###


# In[112]:


PCs = 20
res = 0.3
top_genes=2000
dataset = 'R-AG'
h5adpath = f'../0.loom/{dataset}_for_subset.h5ad'
celltypelist = ['Luminal epithelial cells','Basal epithelial cells','Proliferating epithelial cells','Keratinized epithelial cells']
vars_use=['sample']
NewCellType= {
    "Myoepithelial Cells (Contractile)": [0],        # cluster 0：表达 Cst3、Smtn、Lpp、Actg2，典型收缩蛋白，代表成熟肌上皮细胞
    "Basal-like Epithelial Cells": [1],              # cluster 1：Krt17 等基因提示其为基底上皮细胞，但不具备明显干性
    "Luminal Progenitor/Stem Cells": [2],            # cluster 2：代谢及线粒体相关基因较高，符合早期 luminal 干/祖细胞特征
    "Differentiated Luminal Cells": [3],             # cluster 3：Aldoc、Sephs1 等提示分化成熟的分泌性 luminal 细胞
    "Stromal Cells": [4],                            # cluster 4：Naaladl2、Aff3、Pbx3、Nfia，可能代表处于发育调控阶段的间充质细胞
    "Basal/Myoepithelial Stem Cells": [5],           # cluster 5：Krt5、Fabp5 等经典基因，提示其具有基底干细胞性质
    "Macrophages": [6],                              # cluster 6：Lgals3、Sat1、B2m、S100a10，符合固有免疫细胞特征
    "Cycling/Proliferating Cells": [7],              # cluster 7：Ccnb2、Ccnb1、Cdk1 等细胞周期基因高表达，代表活跃增殖细胞
    "Collagen-producing Fibroblasts": [8],           # cluster 8：Col1a1/Col1a2、Dcn 高表达，典型基质成纤维细胞
    "Vascular/Endothelial Cells": [9],               # cluster 9：Gsn、Timp3 等提示其为血管内皮及相关细胞
    "B Cells": [10],                                 # cluster 10：Ighm 等免疫标记，明确属于 B 细胞群
    "Myofibroblasts": [11],                          # cluster 11：含有 Col1a2、Acot13 等，可能兼具基质和部分收缩功能
    "Schwann Cells": [12],                           # cluster 12：Stard13、S100b、Cdh19 表达，符合周围神经胶质细胞特征
    "Muscle Cells": [13],                            # cluster 13：Cadm2、Lama2、Mybpc1 等肌肉相关基因表达
    "Adipocytes": [14]                               # cluster 14：BC048679、Cyp4b1、Lactb2 等提示脂肪细胞特征
}

update_cell_type = True
clusterlist = ['']
Featuregenes = ['Pip','Esr1','Epcam','Lalba','Top2a','Pgr','Prlr','Acta2','Elf5','Tcf4','Krt1','Ar','Pigr','Cd69','Adipoq','Lum','Vim','Ptprc','Lef1','Tpm2','Krt23','Krt10','Faah','Tm4sf1','Ppl','Wnt11','Krtdap','Sbsn','Dsp','Rab25','Aqp3','Shh','Atp1a1','Atp1b1','Procr']
Featuregenes = [
    'Cst3', 'Smtn', 'Lpp', 'Actg2',         # Myoepithelial Cells (Contractile) (cluster 0)
    'Krt17', 'Plcb1', 'Pde4d', 'Spp1',         # Basal-like Epithelial Cells (cluster 1)
    'Sod1', 'Micos10', 'Echs1', 'Cox14',        # Luminal Progenitor/Stem Cells (cluster 2)
    'Aldoc', 'Sephs1', 'Scp2', 'Eci1',           # Differentiated Luminal Cells (cluster 3)
    'Naaladl2', 'Aff3', 'Pbx3', 'Nfia',          # Stromal Cells (cluster 4)
    'Krt5', 'Fabp5', 'Lgals7', 'Klf5',           # Basal/Myoepithelial Stem Cells (cluster 5)
    'Lgals3', 'Sat1', 'B2m', 'S100a10',          # Macrophages (cluster 6)
    'Ccnb2', 'Ccnb1', 'Secisbp2', 'Cdk1',         # Cycling/Proliferating Cells (cluster 7)
    'Col3a1', 'Col1a2', 'Dcn', 'Col1a1',          # Collagen-producing Fibroblasts (cluster 8)
    'Gsn', 'Alkbh3', 'Timp3', 'Pde3a',           # Vascular/Endothelial Cells (cluster 9)
    'Ighm', 'Mef2c', 'Arhgap15', 'Ptprc',         # B Cells (cluster 10)
    'Col1a2', 'Acot13', 'Col3a1', 'Aldoc',        # Myofibroblasts (cluster 11)
    'Stard13', 'Slc35f1', 'S100b', 'Cdh19',       # Schwann Cells (cluster 12)
    'Cadm2', 'Lama2', 'Mybpc1', 'Aopep',          # Muscle Cells (cluster 13)
    'BC048679', 'Acot13', 'Cyp4b1', 'Lactb2'       # Adipocytes (cluster 14)
,'Krt7','Pip','Pigr','Cdh1','Cldn3'
]

random.seed(2024)
np.random.seed(2024)
cells = 5000
doFeatureplot = True
runharmony = True
subset_celltype = True
do_subset_cluster = False
random_subset_cells = False
use_scvi=False


# In[113]:


###################function###############


# In[114]:


def update_cell_annotations(adata, clustering_column, new_cell_type_dict, update_cell_type):
    if not update_cell_type:
        print("细胞注释更新被禁用")
        return adata


    cluster_to_celltype = {}

    for cell_type, clusters in new_cell_type_dict.items():
        for cluster in clusters:
            cluster_to_celltype[str(cluster)] = cell_type
    adata.obs['subtype'] = adata.obs[clustering_column].map(cluster_to_celltype)

    print("finish updata_cell_annotation")
    print(adata.obs['subtype'].value_counts())
    print(adata.obs['celltype'].value_counts())
    return adata


# In[115]:


def run_harmony(adata,vars_use=vars_use):
    print('running harmony')
    pca_result = adata.obsm['X_pca']
    ho = hm.run_harmony(pca_result, adata.obs, vars_use,random_state=42)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print('finished harmony')
    return adata


# In[116]:


def random_subset(adata,n_cells):
    print(f"random subset cells {n_cells}")
    random_indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)
    adata_sub = adata[random_indices, :]
    return adata_sub


# In[117]:


def export_deg_result(adata):
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


# In[118]:


def run_scvi(adata,batch='stage',res=res):
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    model.train()
    SCVI_LATENT_KEY = "X_scVI"
    adata.obsm[SCVI_LATENT_KEY] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep=SCVI_LATENT_KEY)
    sc.tl.leiden(adata,resolution = res)
    SCVI_MDE_KEY = "X_scVI_MDE"
    adata.obsm[SCVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCVI_LATENT_KEY], accelerator="cpu")
    return(adata)


# In[119]:


def run_preprocess(adata,top_genes):
    # Normalizing to median total counts
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)
    adata.raw=adata
    adata.layers["normalized"] = adata.X.copy()
    print('Finish normalized')
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
    sc.pl.highly_variable_genes(adata)
    print('Finish Varible genes')


# In[120]:


def run_reduceDimension(ea,use_scvi,runharmony,PCs,res):
    import scanpy as  sc
    sc.pp.scale(ea,max_value=10)
    sc.tl.pca(ea,mask_var="highly_variable")
    sc.pl.pca_variance_ratio(ea,log = False)
    if use_scvi:   
        import os
        import tempfile
        import scanpy as sc
        import scvi
        import seaborn as sns
        import torch
        from scib_metrics.benchmark import Benchmarker
        scvi.settings.seed = 0
        print("Last run with scvi-tools version:", scvi.__version__)
        ea=run_scvi(ea)
        sc.pl.embedding(
            ea,
            basis=SCVI_MDE_KEY,
            color=['stage', "leiden",'celltype'],
            frameon=False,
            ncols=1,
            save=f'{dataset}_SCVI.png',
        )
    else:
        print('skip scvi')
    if runharmony:
        run_harmony(ea)
        sc.pp.neighbors(ea, use_rep='X_pca_harmony',n_pcs=PCs,random_state=42,n_neighbors=15)
        print('Finish neighbors')
    else:
        sc.pp.neighbors(ea,n_pcs=PCs,random_state=42)
        print('Finish neighbors')
    sc.tl.leiden(ea,resolution = res,random_state=42)
    print('Finish clustering')
    sc.tl.umap(ea,random_state=42)
    print('Finish UMAP')


# In[121]:


def do_umap_plots(ea,dataset,Featuregenes,doFeatureplot):
    sc.settings.figdir = ''
    sc.pl.umap(
        ea,
        color=["subtype", "celltype", "stage",'sample','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_cluster.png',
    color_map='viridis'
    )

    if doFeatureplot:
        sc.pl.umap(
            ea,
            color=[*Featuregenes,'celltype',"stage",'sample','leiden','gland'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
        ncols=4,
        save=f'{dataset}_cluster_feature.png',
        color_map='plasma'
        )
        print('Featureplot Finished')
    else:
        print('666666')


# In[122]:


def do_DEG(ea,dataset):
    #### Different gene test###
    sc.tl.rank_genes_groups(ea,groupby='leiden',method = 'wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ea,groupby='leiden',n_genes=5,save=f'{dataset}_dotplot_leiden.png',min_logfoldchange=0.25)
    df1 = export_deg_result(adata=ea)
    df1.to_csv(f'{dataset}_ranked_genes_leiden.csv', index=False)


# In[123]:


def do_cell_barplot(ea,dataset):
    celltype_counts = ea.obs.groupby(['stage', 'newcelltype']).size().unstack(fill_value=0)
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
    
    ax.set_title('Percentage of Each Celltype Over Time', fontsize=16)
    ax.set_xlabel('Stage', fontsize=14)
    ax.set_ylabel('Percentage', fontsize=14)
    
    ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.show()
    fig3.savefig(f"{dataset}_raw_bar_plot.png")


# In[124]:


###########################


# In[95]:


#h5adpath=f"../{dataset}.h5ad"


# In[125]:


#h5adpath=f"D:/111/{dataset}_for_subset.h5ad"


# In[126]:


#h5adpath=f"D:/111/{dataset}.h5ad"
h5adpath


# In[127]:


ea_raw=sc.read_h5ad(h5adpath)


# In[128]:


ea_raw.obs['newcelltype'] = ea_raw.obs['celltype']


# In[128]:


ea_raw.obs['subtype'] = ea_raw.obs['celltype']


# In[129]:


if subset_celltype:
    ea = ea_raw[ea_raw.obs['celltype'].isin(celltypelist)].copy()
else:
    ea = ea_raw


# In[130]:


if random_subset_cells:
    print(f'begin subset cells, cells number: {cells}')
    ea = random_subset(ea,n_cells = cells)
else:
    print('666')


# In[131]:


ea


# In[132]:


run_preprocess(ea,top_genes=top_genes)


# In[133]:


run_reduceDimension(ea,use_scvi=use_scvi,runharmony=runharmony,PCs=PCs,res=res)


# In[134]:


if do_subset_cluster:
    ea = ea[~ea.obs['leiden'].isin(clusterlist)].copy()


# In[135]:


ea = update_cell_annotations(ea, clustering_column='leiden', new_cell_type_dict=NewCellType, update_cell_type=update_cell_type)


# In[28]:


ea.write(f'{dataset}_for_DEG.h5ad')


# In[136]:


do_umap_plots(ea=ea,dataset=dataset,Featuregenes=Featuregenes,doFeatureplot=doFeatureplot)


# In[ ]:


do_DEG(ea=ea,dataset=dataset)


# In[143]:


do_cell_barplot(ea=ea,dataset=dataset)


# In[ ]:




