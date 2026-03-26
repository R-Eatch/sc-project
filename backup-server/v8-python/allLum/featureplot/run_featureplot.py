#!/usr/bin/env python
# coding: utf-8

# In[22]:


import scanpy as sc
import gc
sc.set_figure_params(
    dpi=100,          # 笔记本/预览分辨率
    dpi_save=300,     # 导出分辨率，期刊常用 300–600
    vector_friendly=True,  # 让散点在 PDF/SVG 中栅格化（默认就是 True）
    format='pdf'
)
##############################
genes_to_plot =["Epcam","Esr1","Elf5","Lef1","Krt5","Stc2","Atp1b1","Krt7",'Elf3','Runx2',"Gata3","Pdk4",'Pip','Clu']
genes_to_plot =[
    "Lpl",
    "Cldn3",
    "Btn1a1",
    "S100a1",
    "Spp1",
    "Xbp1",
    "Cldn8"
]
obs_list=['newcelltype','stage']
obs_list=[]
##############################
datasetlist = [ 'M-MG',
                'R-MG',
                'S-MG', 
                'R-AG',
                'R-CG',
                'S-AG'
              ]
data_path = "/data01/sunxuebo/project/scrnaseq/v8-python/"


# In[ ]:


for dataset in datasetlist:
    print(f"正在处理数据集：{dataset}")
    adata = sc.read_h5ad(f"{data_path}{dataset}/1.subset/{dataset}_cleaned.h5ad")
    adata.X=adata.layers['normalized']
    valid_genes = [gene for gene in genes_to_plot if gene in adata.var_names or gene == 'newcelltype' or gene == 'stage']
    valid_genes_for_umap=valid_genes+obs_list
    sc.pl.umap(
        adata,
        color=valid_genes_for_umap,
        color_map='viridis',
        ncols=5,
        frameon=False,
        #legend_loc='on data',
        save=f"_{dataset}_featureplot.pdf",
        show=False
    )
    sc.pl.umap(
        adata,
        color=valid_genes_for_umap,
        color_map='viridis',
        ncols=2,
        frameon=False,
        #legend_loc='on data',
        save=f"_{dataset}_featureplot.png",
        show=False
    )
    sc.pl.heatmap(adata, valid_genes_for_umap, groupby='newcelltype', swap_axes=True,save=f"_{dataset}_heatmap.png")
    # sc.pl.embedding(
    # adata,
    # basis='umap',         # Specify which embedding to use (essential!)
    # color='pseudotime',
    # color_map='plasma',
    # frameon=False,
    # colorbar_loc=None,
    # title=f'{dataset}_pseudotime',# This is the key parameter to remove the colorbar
    # save=f"_{dataset}_pseudotime-umap.pdf", # Consider a new filename
    # show=True
    # # You might experiment with the 'size' parameter if points are too large/small
    # # size=10 # Example: Adjust marker size (default depends on data size)
# )
    for obs1 in obs_list:
        sc.pl.violin(
        adata,
        keys= valid_genes,
        groupby=obs1, 
        jitter=0.4,
        rotation=90,
        size=2,stripplot=False,save=f"_{dataset}_{obs1}_violin.png"
        )
    fig=sc.pl.stacked_violin(adata,var_names=valid_genes, groupby='newcelltype', dendrogram=False,return_fig=True,show=False)
    #fig.subtitle(f"{dataset} stacked violin", fontsize=16)
    fig.savefig(f"figures/{dataset}_violin_all.png",dpi=300,bbox_inches="tight")
    fig.savefig(f"figures/{dataset}_violin_all.pdf",dpi=300,bbox_inches="tight")
    del adata
    gc.collect()
print("全部绘制完成，带标签的图片已保存至figures目录中。")

