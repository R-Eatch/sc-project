#!/usr/bin/env python
# coding: utf-8

# In[6]:


import os
os.environ['PYTHONHASHSEED'] = '0'
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
dataset = 'S-MAG'
h5adpath = f'../1.harmony/{dataset}_cleaned.h5ad'
celltypelist = ['ASC','MaSC']
vars_use=['sample']
NewCellType = {
    "StemCells": [2,4,7,6,10]
}
update_cell_type = False
clusterlist = ['']
Featuregenes = ['Lef1','Krt5','Prlr','Fabp5']
random.seed(2024)
np.random.seed(2024)
cells = 5000
doFeatureplot = True
runharmony = False
subset_celltype = True
do_subset_cluster = False
random_subset_cells = False
use_scvi=False


# In[113]:


###################function###############


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
        color=["newcelltype", "celltype", "stage",'sample','leiden','gland'],
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
            color=[*Featuregenes,'newcelltype',"stage",'sample','leiden','gland'],
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
    sc.tl.rank_genes_groups(ea,groupby='newcelltype',method = 'wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ea,groupby='newcelltype',n_genes=10,save=f'{dataset}_dotplot_leiden.png',min_logfoldchange=0.25)
    df1 = export_deg_result(adata=ea)
    df1.to_csv(f'{dataset}_ranked_genes_newcelltype.csv', index=False)


# In[4]:


def do_common_genes(adata, groupby, pval_cutoff=0.05, lfc_cutoff=1.0, out_csv="common_genes.csv"):

    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method='wilcoxon',
        reference='rest'
    )

    # 2. 获取所有指定群体的差异表达结果
    #    当group=None时，会返回所有做了rank_genes_groups的群体的整合结果
    df_all = sc.get.rank_genes_groups_df(adata, group=None)

    # 3. 分别在两个群（例如 group1 与 group2）中，筛选显著上调基因
    #    （你也可以改成循环方式，对多个群取交集）
    group1, group2 = celltypelist
    df_sig_group1 = df_all[
        (df_all['group'] == group1) &
        (df_all['pvals_adj'] < pval_cutoff) &
        (df_all['logfoldchanges'] > lfc_cutoff)
    ]
    df_sig_group2 = df_all[
        (df_all['group'] == group2) &
        (df_all['pvals_adj'] < pval_cutoff) &
        (df_all['logfoldchanges'] > lfc_cutoff)
    ]

    # 4. 求交集
    common_genes = set(df_sig_group1['names']).intersection(df_sig_group2['names'])

    # 5. 筛选原始df_all中这些“共同富集”的基因行
    df_filtered = df_all[df_all['names'].isin(common_genes)]

    # 6. 保存到CSV文件
    df_filtered.to_csv(out_csv, index=False, encoding='utf-8')

    print(f"Common significantly expressed genes (p<{pval_cutoff}, logFC>{lfc_cutoff}):")
    print(common_genes)

    # 若需要进一步使用，函数也可返回 df_filtered
    return df_filtered


# In[5]:


def compare_two_groups_dotplot(
    adata, 
    groupby, 
    group1, 
    group2, 
    pval_cutoff=0.05, 
    lfc_cutoff=1.0, 
    top_n=20, 
    out_csv="deg_two_groups.csv"
):
    # 1. rank_genes_groups：指定 group1 与 group2 做严格两两比较
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        groups=[group1],    # 只对 group1 做测试
        reference=group2,   # 直接对比 group2
        method='wilcoxon'
    )

    # 2. 从adata中提取 group1 vs. group2 的差异表达结果
    #    注意：sc.get.rank_genes_groups_df(adata, group=group1) 会返回group1相对于group2的DEG结果
    df = sc.get.rank_genes_groups_df(adata, group=group1)
    
    # 3. 筛选出显著且满足logFC阈值的基因
    df_filtered = df[
        (df['pvals_adj'] < pval_cutoff) &
        (df['logfoldchanges'] > lfc_cutoff)
    ].copy()

    # 4. 根据logFC从大到小排序，并取前 top_n 个基因用于可视化
    df_filtered.sort_values('logfoldchanges', ascending=False, inplace=True)
    top_genes = df_filtered['names'].head(top_n).tolist()

    # 5. 保存所有筛选后的基因结果到 CSV 文件
    df_filtered.to_csv(out_csv, index=False, encoding='utf-8')
    print(f"保存完成：共有 {df_filtered.shape[0]} 个基因满足阈值，已保存到 {out_csv}")

    # 6. 绘制DotPlot（只针对 group1、group2 这两个类别）
    #    注意：若想看所有群体中这些基因的表达情况，可直接在原始adata上作图
    adata_sub = adata[adata.obs[groupby].isin([group1, group2])]  # 仅保留这两个群
    sc.pl.dotplot(
        adata_sub, 
        var_names=top_genes, 
        groupby=groupby,save='C1vsC2.png'
    )

    # 如果后面还想进一步操作，可以return相关信息
    return df_filtered, top_genes


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
ea_raw.X=ea_raw.layers['counts']


# In[129]:


if subset_celltype:
    ea = ea_raw[ea_raw.obs['newcelltype'].isin(celltypelist)].copy()
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


# In[28]:


#ea.write(f'{dataset}_for_DEG.h5ad')


# In[136]:


do_umap_plots(ea=ea,dataset=dataset,Featuregenes=Featuregenes,doFeatureplot=doFeatureplot)


# In[ ]:


#do_common_genes(adata=ea,groupby='newcelltype')


# In[ ]:


#compare_two_groups_dotplot(adata=ea,groupby='newcelltype',group1='ASC',group2='MaSC')


# In[ ]:


do_DEG(ea=ea,dataset=dataset)


# In[ ]:




