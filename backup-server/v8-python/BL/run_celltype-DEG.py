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
import loompy
print('import successful')


# In[2]:


###global variable###


# In[107]:


dataset = 'R-MG'
h5adpath=f'{dataset}_for_DEG.h5ad'
NewCellType = {
    "StemCells": [1],
    "LumHR": [0,7],
    "Lum-Adipoq": [6],
    "Lum-basal": [1,8],
    "LumSEC":[2,3],
    "Stem-Lum":[5]
}
update_cell_type = False
Featuregenes = ['Esr1','Epcam','Top2a','Acta2','Elf5','Tcf4','Krt1','Prlr','Ar','Pigr','Vim','Anxa1','Krt5','Adipoq','Pgr','Cd74']
genes= Featuregenes
group='newcelltype'
random.seed(2024)
np.random.seed(2024)
random_state = 2024
clusterlist = ['']
celltypelist = ['']
cells = 5000
doFeatureplot = True
subset_celltype = False
do_subset_cluster = False
random_subset_cells = False


# In[4]:


###################function###############


# In[5]:


def update_cell_annotations(adata, clustering_column, new_cell_type_dict, update_cell_type):
    if not update_cell_type:
        print("细胞注释更新被禁用")
        return adata


    cluster_to_celltype = {}

    for cell_type, clusters in new_cell_type_dict.items():
        for cluster in clusters:
            cluster_to_celltype[str(cluster)] = cell_type
    adata.obs['newcelltype'] = adata.obs[clustering_column].map(cluster_to_celltype)

    print("finish updata_cell_annotation")
    print(adata.obs['newcelltype'].value_counts())
    print(adata.obs['celltype'].value_counts())
    return adata


# In[7]:


def random_subset(adata,n_cells):
    print(f"random subset cells {n_cells}")
    random_indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)
    adata_sub = adata[random_indices, :]
    return adata_sub


# In[8]:


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


# In[82]:


def do_umap_plots(ea,dataset,Featuregenes,doFeatureplot):
    sc.settings.figdir = ''
    sc1=sc.pl.umap(
        ea,
        color=["newcelltype", "celltype", "stage",'sample','leiden','gland'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
        return_fig=True,
    color_map='viridis'
    )
    sc1.savefig(f'{dataset}_celltype.pdf')
    sc1.savefig(f'{dataset}_celltype.png')
    sc2=sc.pl.umap(
    ea,
    color=["newcelltype"],
    wspace=0.5,
    size=3,
        return_fig=True
    )
    sc2.savefig(f'{dataset}_newcelltype.pdf')
    sc2.savefig(f'{dataset}_newcelltype.png')
    if doFeatureplot:
        fp=sc.pl.umap(
            ea,
            color=[*Featuregenes,'newcelltype','celltype','stage','sample','leiden','gland','species'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
        ncols=4,
        save=f'{dataset}_feature.png',
        return_fig=True,
        color_map='viridis'
        )
        fp.savefig(f'{dataset}_feature.pdf')
        fp.savefig(f'{dataset}_feature.png')
        
        print('Featureplot Finished')
    else:
        print('666666')


# In[13]:


def do_DEG(ea,dataset):
    #### Different gene test###
    sc.tl.rank_genes_groups(ea,groupby='newcelltype',method = 'wilcoxon')
    sc.pl.rank_genes_groups_dotplot(ea,groupby='newcelltype',n_genes=7,save=f'{dataset}_dotplot.png',min_logfoldchange=0.25)
    df1 = export_deg_result(adata=ea)
    df1.to_csv(f'{dataset}_ranked_genes.csv', index=False)


# In[93]:


def do_marker_plot(adata, genes, group, dataset, color_map="RdBu_r"):

    # 绘制 DotPlot
    dotplot = sc.pl.dotplot(
        adata,
        var_names=genes,
        groupby=group,
        swap_axes=True,         # 纵向显示
        color_map=color_map,    # 配色方案
        dot_max=1.0,            # 可根据需要调整点的最大大小
        standard_scale='var',   # 标准化基因数据
        return_fig=True,
        title=f'Marker Genes Across Different {group}'
        # 返回图对象
    )
    
    # 保存图像为 PNG 和 PDF
    png_path = f"{dataset}_dotplot_marker.png"
    pdf_path = f"{dataset}.dotplot_marker.pdf"
    dotplot.savefig(png_path)
    dotplot.savefig(pdf_path)
    
    print(f"DotPlot saved as {png_path} and {pdf_path}")
    # 显示图表
    plt.show()
    sc.pl.violin(
    adata,
    keys=genes,
    groupby=group,  # 替换为分组列名
    jitter=0.4,
    scale="width",
    rotation=45,
    size=2,
    stripplot=False,save=f"{dataset}.violin_marker.pdf")
    sc.pl.violin(
    adata,
    keys=genes,
    groupby=group,  
    jitter=0.4,
    scale="width",
    rotation=45,
    size=2,
    stripplot=False,save=f"{dataset}.violin_marker.png"
    )


# In[127]:


def do_cell_barplot(ea, dataset):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns
    
    # 计算细胞类型比例
    celltype_counts = ea.obs.groupby(['stage', 'newcelltype']).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    celltype_percentages_long = celltype_percentages.reset_index().melt(id_vars='stage', var_name='Celltype', value_name='Percentage')
    
    # 默认配色
    default_palette = ea.uns['newcelltype_colors']
    
    # 创建绘图
    times = celltype_percentages.index
    cell_types = celltype_percentages.columns
    colors = default_palette

    bottom = pd.Series([0] * len(times), index=times)

    fig, ax = plt.subplots(figsize=(16, 8))
    
    # 设置背景为白色并去除网格
    ax.set_facecolor('white')
    sns.set(style="white")
    
    # 绘制柱状图，调整柱子的宽度
    bar_width = 0.25  # 调整柱子宽度
    # 生成 X 轴的位置，增大间距
    x_positions = np.arange(len(times)) * 0.5  # 1.5 为间距比例，可根据需求调整
    
    for cell_type, color in zip(cell_types, colors):
        percentages = celltype_percentages[cell_type]
        ax.bar(x_positions, percentages, bottom=bottom, color=color, label=cell_type, width=bar_width)
        bottom += percentages
    
    # 设置 X 轴标签的位置
    ax.set_xticks(x_positions)
    ax.set_xticklabels(times, rotation=45)

    # 设置标题和坐标轴标签
    ax.set_title(f'Percentage of Each Celltype Over Stage in {dataset}', fontsize=16)  # 修改标题的"TIME"为"stage"
    ax.set_xlabel('Stage', fontsize=14)
    ax.set_ylabel('Percentage', fontsize=14)

    # 设置图例
    ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # 调整x轴标签的角度
    plt.xticks(rotation=45)
    
    # 布局调整
    plt.tight_layout()
    
    # 显示图表
    plt.show()
    
    # 保存图表
    fig.savefig(f"{dataset}_bar_plot.png")
    fig.savefig(f"{dataset}_bar_plot.pdf")


# In[15]:


###########################


# In[16]:


h5adpath=f'{dataset}_for_DEG.h5ad'


# In[108]:


#h5adpath=f"D:/111/{dataset}_cleaned.h5ad"


# In[98]:


#h5adpath=f"D:/111/{dataset}.h5ad"
h5adpath


# In[109]:


ea = ad.read_h5ad(h5adpath)


# In[28]:


if do_subset_cluster:
    ea = ea[~ea.obs['leiden'].isin(clusterlist)].copy()


# In[29]:


ea = update_cell_annotations(ea, clustering_column='leiden', new_cell_type_dict=NewCellType, update_cell_type=update_cell_type)


# In[ ]:


if subset_celltype:
    ea = ea[~ea.obs['newcelltype'].isin(celltypelist)].copy()


# In[30]:


ea.write(f'{dataset}_cleaned.h5ad')


# In[111]:


do_umap_plots(ea=ea,dataset=dataset,Featuregenes=Featuregenes,doFeatureplot=doFeatureplot)


# In[32]:


do_DEG(ea=ea,dataset=dataset)


# In[128]:


do_cell_barplot(ea=ea,dataset=dataset)


# In[106]:


do_marker_plot(adata=ea,genes=genes,group=group,dataset=dataset)


# In[ ]:




