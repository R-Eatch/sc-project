#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as  sc
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import random
import seaborn as sns
import scanpy.external as sce
import harmonypy as hm
print('import successful')
sc.settings.n_jobs = 1


# In[4]:


###global variable###
PCs = 20 
res = 0.2
dataset='ALL-MG'
if dataset == "ALL-MG":
    dataset1='M-MG'
    dataset2='R-MG'
    dataset3='S-MG'
    NewCellType = {
    "Basal":["Basal"],
    "Luminal HR":["LumHR"],
    "Luminal SEC":["LumSEC-Lac", "LumSEC-Lip", "LumSEC-Mgp", "LumSEC-Vcam1"],
    "Orther Epithelial":["Epi-Fibro", "Epi-Krt7", "Epi-Lgals7","Lum-Basal", "Lum-Immune"],
    "MaSC":["MaSC", "MaSC-Pro", "MaSC-t2-sg", "Lum-Kit"]
    }
    Featuregenes = ['Krt5','Acta2','Krt17','Prlr','Pgr','Esr1','Epcam','Elf5','Lalba','Lef1','Krt5','Gli2']
elif dataset == 'ALL-AG':
    dataset1='R-AG'
    dataset2='R-CG'
    dataset3='S-AG'
    NewCellType = {
    "Basal":["Basal"],
    "Luminal SEC":["LumSEC-AG-Pip", "LumSEC-AG-t1", "LumSEC-AG-t2", "LumSEC-Lip-CG","Lum-Stat4", "Lum-Tm4sf4"],
    "Orther Epithelial":["Epi-Fibro", "Epi-Hmgcr", "Epi-Krt7", "Epi-Lalba", "Epi-Lgals7",
    "Epi-Lip", "Epi-Pro"],
    "ASCs":["ASC-rb", "ASC-sg", "CSC"]
    }
    Featuregenes = ['Krt5','Acta2','Krt17','Epcam','Pip','Sephs1','Krt23','Mif','Melk']
path1=f"../../{dataset1}/1.subset/{dataset1}_cleaned.h5ad"
path2=f"../../{dataset2}/1.subset/{dataset2}_cleaned.h5ad"
path3=f"../../{dataset3}/1.subset/{dataset3}_cleaned.h5ad"
vars_use = ['batchlist']
random.seed(12345)
np.random.seed(12345)

doFeatureplot = True
runharmony = True
do_subset_cluster = False
subset_celltype = False
re_find_hvgenes = True
update_cell_type = True
random_state=2024
celltypelist = ['Lum-Immun-like','Lum-Ker-like','Lum-Fibro-like']
sc.settings.figdir = ''
#########################


# In[5]:


####function######


# In[ ]:


def run_harmony(adata,vars_use=vars_use):
    print('running harmony')
    pca_result = adata.obsm['X_pca']
    ho = hm.run_harmony(pca_result, adata.obs, vars_use,random_state=42)
    adata.obsm['X_pca_harmony'] = ho.Z_corr.T
    print('finished harmony')
    return adata


# In[7]:


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


# In[8]:


def subsetcelltype(adata,celltypelist):
    adata = adata[~adata.obs['newcelltype'].isin(celltypelist),:].copy()
    return adata


# In[ ]:


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

    adata.obs['mergedcelltype'] = adata.obs[clustering_column].map(cluster_to_celltype)
    print("Finished updating cell annotation:")
    print(adata.obs['mergedcelltype'].value_counts())
    return adata


# In[ ]:


#######################


# In[ ]:


adata1 = sc.read_h5ad(path1)
adata2 = sc.read_h5ad(path2)
adata3 = sc.read_h5ad(path3)


# In[ ]:


adata1.obs = adata1.obs.drop(columns=["cellid"])
adata2.obs = adata2.obs.drop(columns=["cellid"])
adata3.obs = adata3.obs.drop(columns=["cellid"])


# In[ ]:


if subset_celltype:
    adata1 = subsetcelltype(adata1,celltypelist)
    adata2 = subsetcelltype(adata2,celltypelist)
    unique_gland = '-'.join(adata1.obs['gland'].unique().tolist())
    adata1.obs['newcelltype'] = [unique_gland + '-' + str(newtype) for newtype in adata1.obs['newcelltype']]
    unique_gland = '-'.join(adata2.obs['gland'].unique().tolist())
    adata2.obs['newcelltype'] = [unique_gland + '-' + str(newtype) for newtype in adata2.obs['newcelltype']]


# In[11]:


adata_combined = sc.concat(
    [adata1, adata2, adata3],
    join='outer',              # 保留所有变量（默认 inner）
    label='batchlist',        # 新增 obs['glandbatch'] 列
    keys=['batch1','batch2','batch3']  # 每个对象对应的类别标签
)


# In[ ]:


adata_combined.obs['mergedcelltype']=adata_combined.obs['newcelltype']
adata_combined=update_cell_annotations(adata_combined,clustering_column='newcelltype',new_cell_type_dict=NewCellType,update_cell_type=update_cell_type)


# In[13]:


if re_find_hvgenes:
    adata_combined.X = adata_combined.layers['normalized'].copy()
    # sc.pp.normalize_total(adata_combined)
    # sc.pp.log1p(adata_combined)
    # print('Finish normalized')
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=2000)
    print('Finish Varible genes')
    sc.pp.scale(adata_combined,max_value=10)


# In[14]:


sc.tl.pca(adata_combined,use_highly_variable=True,random_state=42)
sc.pl.pca_variance_ratio(adata_combined,log = False)


# In[27]:


if runharmony:
    run_harmony(adata_combined)
    sc.pp.neighbors(adata_combined, use_rep='X_pca_harmony',n_pcs=PCs,random_state=random_state,n_neighbors=15)
else:
    sc.pp.neighbors(adata_combined,n_pcs= PCs,random_state=random_state)
#sc.tl.louvain(ea,resolution = res)
sc.tl.leiden(adata_combined,resolution = res,random_state=random_state)
print('Finish clustering')
sc.tl.umap(adata_combined,random_state=random_statem)
print('Finish UMAP')


# In[ ]:





# In[29]:


sc.pl.umap(
        adata_combined,
        color=["newcelltype", "celltype", "mergedcelltype",'stage','sample','species'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_merged_celltype.png',
    color_map='viridis'
    )


# In[29]:


sc.pl.umap(
        adata_combined,
        color=["newcelltype", "celltype", "mergedcelltype",'stage','sample','species'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_merged_celltype.pdf',
    color_map='viridis'
    )


# In[ ]:


adata_combined.write(f'{dataset}_cleaned.h5ad')


# In[30]:


if doFeatureplot:
    sc.pl.umap(
        adata_combined,
        color=[*Featuregenes,'newcelltype','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save=f'{dataset}_feature.pdf',
    color_map='seismic',
    )
    print('Featureplot Finished')
else:
    print('666666')


# In[ ]:


#### Different gene test###
sc.tl.rank_genes_groups(adata_combined,groupby='mergedcelltype',method = 'wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata_combined,groupby='mergedcelltype',n_genes=8,save=f'{dataset}_dotplot.png',min_logfoldchange=0.1)


# In[ ]:


df1 = export_deg_result(adata=adata_combined)


# In[ ]:


df1.to_csv(f'{dataset}_celltype_genes.csv', index=False)


# In[ ]:


sc.tl.rank_genes_groups(adata_combined,groupby='stage',method = 'wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata_combined,groupby='stage',n_genes=8,save=f'{dataset}_dotplot_leiden.png',min_logfoldchange=0.1)


# In[ ]:


df1 = export_deg_result(adata=adata_combined)
df1.to_csv(f'{dataset}_stage_genes.csv', index=False)


# In[ ]:


ea = adata_combined


# In[ ]:


ea.obs['stage'] = ea.obs['stage'].astype(str)
celltype_counts = ea.obs.groupby(['stage', 'mergedcelltype']).size().unstack(fill_value=0)
celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
celltype_percentages_long = celltype_percentages.reset_index().melt(id_vars='stage', var_name='Celltype', value_name='Percentage')


# In[ ]:


default_palette = ea.uns['mergedcelltype_colors']
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
ax.set_xlabel('stage', fontsize=14)
ax.set_ylabel('Percentage', fontsize=14)

ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xticks(rotation=45)
plt.tight_layout()
plt.show()


# In[ ]:


fig3.savefig(f"{dataset}_bar_plot.png")

