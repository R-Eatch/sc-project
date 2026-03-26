#!/usr/bin/env python
# coding: utf-8

# In[4]:


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


# In[5]:


###global variable###


# In[74]:


dataset = 'S-EG'
h5adpath=f'{dataset}_cleaned.h5ad'
NewCellType = {
    "Basal": [0],
}
use_spefic_color= False
update_cell_type = False
group='newcelltype'
Featuregenes = ['Esr1','Epcam','Top2a','Acta2','Elf5','Tcf4','Krt1','Prlr','Ar','Pigr','Vim','Anxa1','Krt5','Adipoq','Pgr','Cd74','Lalba','Dcn','Col1a1']
genes= [
        'Prlr','Krt23','Top2a',#ASC
        #'Barx2','Meis2','Klk6',#CSC
        #'Acta2','Tpm2','Krt5', #basal
        #'Adipoq','Nr1h3','Gpat3', #adi
        #'Col1a1','Lum','Vim', #fibro
        #'Krt1','Krt10','Krtdap', #krt
        #'Cd69','Ptprc','Cd86',#IMMUNE
        #'Esr1','Pgr','Prlr', # MG-HR
        #'Lalba','Elf5','Epcam',#MG-sec
        #'Esr1','Pip','Pigr', # AG-Esr+
        #'Faah','Slc27a4','Hacd2',#AG-sec1-lip
        #'Tm4sf1','Wnt11','Ppl', #AG-SEC2-ker
        #'Top2a','Lef1','Meis2',#Masc
        #'Lalba','Elf5','Epcam',#AG-SEC-MG-LIKE
        #'Acot7','Elovl6','Fabp3',#CG-SEC-lip
        #'Tnfrsf21','Erbb4','Creb5',#CG-reg
       ]
random.seed(2024)
np.random.seed(2024)
random_state = 2024
clusterlist = ['']
celltypelist = ['']
cells = 5000
doFeatureplot = False
subset_celltype = False
do_subset_cluster = False
random_subset_cells = False

color_dict = {
    "MaSC": "#a6cee3",
    "Lum-Ker-like": "#1f78b4",
    "Lum-Imm-like": "#ffd700",
    "Lum-Fib-like": "#000000",
    "LumHR-MG": "#33a02c",
    "LumSEC-MG": "#e31a1c",
    "LumSEC-AG-t1": "#fdbf6f",
    "LumSEC-AG-t2": "#6a3d9a",
    "ASC": "#cab2d6",
    "CSC": "#ff7f00",
    "LumSEC-MG-like": "#fb9a99",
    "Lum-Esr+":"#b15928",
    "LumSEC-CG": "#b2df8a",
    "Lum-t2": "#808080",
    "Basal": "#0047AB",
    "Lum-Adi-like": "#FF5733"
}


# In[7]:


###################function###############


# In[8]:


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


# In[9]:


def random_subset(adata,n_cells):
    print(f"random subset cells {n_cells}")
    random_indices = np.random.choice(adata.n_obs, size=n_cells, replace=False)
    adata_sub = adata[random_indices, :]
    return adata_sub


# In[10]:


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


# In[11]:


def do_celltype_color(ea, color_dict):
    unique_cell_types = ea.obs['newcelltype'].unique()
    filtered_color_dict = {k: color_dict[k] for k in unique_cell_types if k in color_dict}
    sorted_filtered_color_dict = {
        k: filtered_color_dict[k]
        for k in sorted(filtered_color_dict.keys())  # 对过滤后的键排序
    }
    print(f"Filtered color dictionary (sorted): {sorted_filtered_color_dict}")
    filtered_palette = list(sorted_filtered_color_dict.values())
    sc.pl.umap(ea, color='newcelltype', palette=filtered_palette)


# In[12]:


def do_umap_plots(ea,dataset,Featuregenes,genes,doFeatureplot):
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
        ncols=3,
        save=f'{dataset}_feature.png',
        return_fig=True,
        color_map='viridis'
        )
        fp.savefig(f'{dataset}_feature.pdf')
        fp.savefig(f'{dataset}_feature.png')
    if doFeatureplot:
        fp=sc.pl.umap(
            ea,
            color=[*genes,'newcelltype','celltype','stage','sample','leiden','gland','species'],
            # increase horizontal space between panels
            wspace=0.5,
            size=3,
        ncols=3,
        return_fig=True,
        color_map='viridis'
        )
        fp.savefig(f'{dataset}_markergenes.pdf')
        fp.savefig(f'{dataset}_markergenes.png')
        
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


# In[14]:


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


# In[15]:


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


# In[16]:


def do_cell_barplot(ea, dataset,celltype_name):
    import matplotlib.pyplot as plt
    import pandas as pd
    import seaborn as sns

    # 计算细胞类型比例
    celltype_counts = ea.obs.groupby(['stage', f'{celltype_name}']).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100
    celltype_percentages_long = celltype_percentages.reset_index().melt(id_vars='stage', var_name='Celltype', value_name='Percentage')

    # 默认配色
    default_palette = ea.uns[f'{celltype_name}_colors']

    # 对 stage 按降序排列
    celltype_percentages = celltype_percentages.sort_index(ascending=False)
    times = celltype_percentages.index  # 更新排序后的 stage 顺序
    cell_types = celltype_percentages.columns
    colors = default_palette

    bottom = pd.Series([0] * len(times), index=times)

    # 设置图像比例为 8:4
    fig, ax = plt.subplots(figsize=(8, 4))
    
    # 绘制横向柱状图，调整柱间距
    bar_width = 0.8  # 设置柱间距为 0
    for cell_type, color in zip(cell_types, colors):
        percentages = celltype_percentages[cell_type]
        ax.barh(times, percentages, left=bottom, color=color, label=cell_type, height=bar_width)  # 横向柱状图
        bottom += percentages

    # 设置标题和坐标轴标签
    ax.set_title(f'Percentage of Each Celltype Over Stage in {dataset}', fontsize=16)
    ax.set_xlabel('Percentage', fontsize=14)
    ax.set_ylabel('Stage', fontsize=14)

    # 设置图例
    ax.legend(title='Celltype', bbox_to_anchor=(1.05, 1), loc='upper left')

    # 调整布局
    plt.tight_layout()

    # 保存图表
    plt.savefig(f"{dataset}_{celltype_name}_bar_plot.png", dpi=300)
    plt.savefig(f"{dataset}_{celltype_name}_bar_plot.pdf")
    plt.show()


# In[17]:


###########################


# In[18]:


#h5adpath=f'{dataset}_for_DEG.h5ad'


# In[351]:


#h5adpath=f"D:/111/{dataset}_cleaned.h5ad"


# In[75]:


#h5adpath=f"D:/111/{dataset}_cleaned.h5ad"
h5adpath


# In[76]:


ea = ad.read_h5ad(h5adpath)


# In[28]:


if do_subset_cluster:
    ea = ea[~ea.obs['leiden'].isin(clusterlist)].copy()


# In[29]:


ea = update_cell_annotations(ea, clustering_column='leiden', new_cell_type_dict=NewCellType, update_cell_type=update_cell_type)


# In[ ]:


if subset_celltype:
    ea = ea[~ea.obs['newcelltype'].isin(celltypelist)].copy()


# In[358]:


if use_spefic_color:
    do_celltype_color(ea,color_dict)


# In[30]:


#ea.write(f'{dataset}_cleaned.h5ad')


# In[260]:


do_umap_plots(ea=ea,dataset=dataset,Featuregenes=Featuregenes,doFeatureplot=doFeatureplot,genes=genes)


# In[32]:


do_DEG(ea=ea,dataset=dataset)


# In[46]:


do_cell_barplot(ea=ea,dataset=dataset,celltype_name='newcelltype')
do_cell_barplot(ea=ea,dataset=dataset,celltype_name='celltype')

# In[77]:


do_marker_plot(adata=ea,genes=genes,group=group,dataset=dataset)


# In[ ]:


##### plot goood fig #######


# In[362]:


import matplotlib.patches as mpatches

def add_color_strip_right_and_below_dotplot(dp, adata, groupby="newcelltype"):
    # 1) 获取主绘图区 & figure
    ax_main = dp.ax_dict["mainplot_ax"]
    fig = ax_main.figure

    # 2) 获取主绘图区在整个 figure 中的坐标 (x0, y0, w, h)
    x0, y0, w, h = ax_main.get_position().bounds

    # 准备在主图右侧留一条 Axes
    cbar_width = 0.03         # 条带宽度（可酌情调节）
    gap = 0.02                # 与主图之间的空隙

    # 右侧位置
    cbar_x = x0 + w + gap

    # 3) 新建一个 Axes，用来画竖着的颜色条带
    cbar_ax_vertical = fig.add_axes([cbar_x, y0, cbar_width, h])  # 右侧条带的坐标
    cbar_ax_vertical.set_xticks([])
    cbar_ax_vertical.set_yticks([])
    cbar_ax_vertical.set_xlim(0, 1)
    cbar_ax_vertical.set_ylim(0, 3)  # 每个颜色长度为 3，按比例调整
    cbar_ax_vertical.axis('off')  # 不显示坐标轴和边框

    # 准备在主图下方留一条 Axes
    cbar_height = 0.03       # 条带高度
    cbar_y = y0 - cbar_height - gap

    cbar_ax_horizontal = fig.add_axes([x0, cbar_y, w, cbar_height])  # 下方条带的坐标
    cbar_ax_horizontal.set_xticks([])
    cbar_ax_horizontal.set_yticks([])
    cbar_ax_horizontal.set_xlim(0, 1)
    cbar_ax_horizontal.set_ylim(0, 1)
    cbar_ax_horizontal.axis('off')  # 不显示坐标轴和边框

    # 确认分类数据和颜色列表
    if not hasattr(adata.obs[groupby], "cat"):
        raise ValueError(f"{groupby} is not a categorical column in adata.obs")
    celltype_order = adata.obs[groupby].cat.categories
    color_list = adata.uns.get(f"{groupby}_colors", None)
    if color_list is None:
        raise ValueError(f"adata.uns does not have key '{groupby}_colors'. Please check your data.")
    if len(color_list) < len(celltype_order):
        raise ValueError("color_list has fewer colors than celltype_order. Check your data.")

    # 竖着的颜色条（顺序反转）
    n_celltypes = len(celltype_order)
    rect_height_vertical = 3.0 / n_celltypes  # 将 cbar_ax_vertical 的 y 方向均分

    for i, ct in enumerate(reversed(celltype_order)):  # 顺序反转
        y_bottom = i * rect_height_vertical
        color_here = color_list[len(celltype_order) - i - 1]
        rect = mpatches.Rectangle((0, y_bottom), 
                                  1,              # 宽度 = 1 (充满 cbar_ax_vertical 的 x 方向)
                                  rect_height_vertical,    # 高度根据 n_celltypes 调整
                                  color=color_here,
                                  linewidth=0
                                  )
        cbar_ax_vertical.add_patch(rect)

    # 下方的颜色条
    rect_width_horizontal = 1.0 / n_celltypes  # 将 cbar_ax_horizontal 的 x 方向均分

    for i, ct in enumerate(celltype_order):
        x_left = i * rect_width_horizontal
        color_here = color_list[i]
        rect = mpatches.Rectangle((x_left, 0), 
                                  rect_width_horizontal, 
                                  1,              # 高度 = 1 (充满 cbar_ax_horizontal 的 y 方向)
                                  color=color_here,
                                  linewidth=0
                                  )
        cbar_ax_horizontal.add_patch(rect)


# In[366]:


if use_spefic_color:
    dp = sc.pl.dotplot(
        ea,
        var_names=genes,
        groupby=group,
        swap_axes=True,         # 纵向显示
        dot_max=1.0,            # 可根据需要调整点的最大大小
        standard_scale='var', return_fig=True,show=False,
        cmap="RdBu_r",colorbar_title="Scaled expression"
        )
    
    ax_main = dp.get_axes()["mainplot_ax"]
    
    ax_main = dp.ax_dict["mainplot_ax"]
    ax_main.set_xticklabels([])
    ax_main.set_xlabel("")
    ax_main.set_title("")
    ax_main.tick_params(axis='x', which='both', bottom=False, top=False)
    #if "color_legend_ax" in dp.ax_dict:
    #    dp.ax_dict["color_legend_ax"].set_visible(False)
    #if "size_legend_ax" in dp.ax_dict:
    #    dp.ax_dict["size_legend_ax"].set_visible(False)
    
    add_color_strip_right_and_below_dotplot(dp, ea, groupby="newcelltype")
    fig1 = ax_main.figure


# In[174]:


def plot_barplot_with_colorbar(ea, dataset, groupby_stage="stage", groupby_celltype="newcelltype"):
    """
    绘制细胞比例 Barplot，并添加从下到上的渐变颜色条（无标注、无刻度）。
    """

    # 1) 计算数据：按 stage 和细胞类型分组，计算百分比
    celltype_counts = ea.obs.groupby([groupby_stage, groupby_celltype]).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100

    # 2) 自定义 stage 排序：确保 stage 是按特定顺序排列
    # 假设 stage 是类似 ["stage1", "stage2", ...] 的字符串，提取数字部分并排序
    stage_order = sorted(celltype_percentages.index, key=lambda x: int(x.replace("stage", "")), reverse=True)
    celltype_percentages = celltype_percentages.reindex(stage_order)

    # 转换 stage 为数值类型，用于渐变色条
    stage_numeric = [int(x.replace("stage", "")) for x in stage_order]

    stages = celltype_percentages.index  # 排序后的 stage
    cell_types = celltype_percentages.columns

    # 获取颜色
    color_list = ea.uns[f"{groupby_celltype}_colors"]  # 假设颜色与 cell_types 一一对应
    if len(color_list) < len(cell_types):
        raise ValueError("颜色列表长度不足，请检查 adata.uns['newcelltype_colors']")

    # 3) 创建一个 Figure 对象
    fig, ax = plt.subplots(figsize=(8, 6))

    # 4) 绘制横向柱状图（调整柱间距）
    bar_width = 0.8  # 柱子的宽度（间距 = 1 - bar_width）
    bottom = np.zeros(len(stages))  # 累积百分比

    for i, cell_type in enumerate(cell_types):
        percentages = celltype_percentages[cell_type]
        color = color_list[i]
        ax.barh(stages, percentages, left=bottom, color=color, height=bar_width, label=cell_type)
        bottom += percentages

    # 5) 图表美化：去掉边框、y 轴刻度、图例等
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)

    ax.set_yticks([])  # 去掉 y 轴刻度
    ax.set_ylabel("")  # 去掉 y 轴标签
    ax.set_xlabel("Percentage (%)", fontsize=12)  # 设置 x 轴标签

    # 如果需要去掉图例
    ax.legend().remove()

    # 6) 添加颜色条表示 stage 的渐变
    # 使用 Normalize 将 stage_numeric 的范围映射到渐变色条
    norm = plt.Normalize(vmin=min(stage_numeric), vmax=max(stage_numeric))
    sm = plt.cm.ScalarMappable(cmap="viridis_r", norm=norm)  # 注意 cmap 用 viridis_r，反转渐变方向
    
    # 在右侧创建一个渐变颜色条
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    cbar = plt.colorbar(sm, cax=cbar_ax, orientation="vertical", extend="min")  # 只在下方延伸为尖状
    
    # 删除颜色条的刻度和标签
    cbar.ax.tick_params(labelsize=0, length=0)  # 隐藏刻度和数字
    cbar.ax.set_title("")  # 删除颜色条标题
    
    # 设置颜色条标签为 "stage"，并垂直显示
    cbar.ax.set_ylabel("stage", rotation=0, labelpad=15)  # 旋转270度并调整间距

    # 7) 返回 Figure 对象
    return fig


# In[175]:


if use_spefic_color:
    fig3 = plot_barplot_with_colorbar(ea, dataset="MyDataset", groupby_stage="stage", groupby_celltype="newcelltype")


# In[170]:


import matplotlib.pyplot as plt
import scanpy as sc

def plot_umap_without_borders(ea, dataset, color="newcelltype"):
    """
    绘制没有边框、标题的 UMAP 图，在左下角添加一个坐标轴。
    返回 Figure 对象。
    """

    # 使用 scanpy 的绘图函数，返回 Figure 对象
    fig = sc.pl.umap(
        ea,
        color=color,
        size=3,
        show=False,          # 不直接展示
        return_fig=True      # 返回 Figure 对象
    )

    # 获取主绘图区
    ax = fig.axes[0]

    # 删除标题
    ax.set_title("")

    # 删除坐标轴刻度和标签
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")

    # 删除边框
    for spine in ax.spines.values():
        spine.set_visible(False)

    # 在左下角添加一个坐标轴
    arrow_length = 0.2  # 坐标轴长度
    arrow_x_start = 0.05  # 箭头起点的 x 坐标（相对轴的范围）
    arrow_y_start = 0.05  # 箭头起点的 y 坐标

    # 绘制箭头（UMAP1）
    ax.annotate(
        "", xy=(arrow_length, 0), xytext=(0, 0),
        arrowprops=dict(facecolor="black", shrink=0, width=1, headwidth=5),
        xycoords="axes fraction", textcoords="axes fraction",
        annotation_clip=False
    )

    # 绘制箭头（UMAP2）
    ax.annotate(
        "", xy=(0, arrow_length), xytext=(0, 0),
        arrowprops=dict(facecolor="black", shrink=0, width=1, headwidth=5),
        xycoords="axes fraction", textcoords="axes fraction",
        annotation_clip=False
    )

    # 添加标签（UMAP1 和 UMAP2）
    ax.text(
        arrow_x_start + 0.2, arrow_y_start, "UMAP1",
        fontsize=10, ha="center", va="center", transform=ax.transAxes
    )
    ax.text(
        arrow_x_start, arrow_y_start + 0.2, "UMAP2",
        fontsize=10, ha="center", va="center", transform=ax.transAxes
    )

    # 返回修改后的 Figure
    return fig


# In[171]:


# 绘制 UMAP 图并返回 Figure 对象
if use_spefic_color:
    fig2 = plot_umap_without_borders(ea, dataset="MyDataset", color="newcelltype")


# In[160]:


import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.image as mpimg
import tempfile
import os

def combine_three_plots(fig1, fig2, fig3, dataset="MyDataset", output_dir="."):
    """
    将 DotPlot (fig1), 细胞类型 UMAP (fig2), 细胞比例 Barplot (fig3) 整合到一张图中。
    - 仅排列 fig1, fig2, fig3
    """
    # 创建一个大图，指定 figsize (可以根据需求调整大小)
    fig = plt.figure(figsize=(20, 14))

    # 使用 GridSpec 定义布局
    gs = GridSpec(2, 2, width_ratios=[2, 1.5], height_ratios=[1.2, 1], wspace=0.3, hspace=0.3)

    # 创建临时目录保存子图
    with tempfile.TemporaryDirectory() as tmpdir:
        # 保存子图为临时图片
        fig1_path = os.path.join(tmpdir, "fig1.png")
        fig2_path = os.path.join(tmpdir, "fig2.png")
        fig3_path = os.path.join(tmpdir, "fig3.png")
        fig1.savefig(fig1_path, dpi=300, bbox_inches="tight")
        fig2.savefig(fig2_path, dpi=300, bbox_inches="tight")
        fig3.savefig(fig3_path, dpi=300, bbox_inches="tight")

        # 加载子图图片并嵌入大图
        # --- 左侧：DotPlot 占用两行 ---
        ax1 = fig.add_subplot(gs[:, 0])
        img1 = mpimg.imread(fig1_path)
        ax1.imshow(img1)
        ax1.axis("off")

        # --- 右上角：细胞类型的 UMAP (放大) ---
        ax2 = fig.add_subplot(gs[0, 1])
        img2 = mpimg.imread(fig2_path)
        ax2.imshow(img2)
        ax2.axis("off")

        # --- 右下角：细胞比例的 Barplot ---
        ax3 = fig.add_subplot(gs[1, 1])
        img3 = mpimg.imread(fig3_path)
        ax3.imshow(img3)
        ax3.axis("off")

        # 设置总标题
        fig.suptitle(f"{dataset}-subcelltype", fontsize=16, y=0.98)

        # 保存整合后的最终图
        combined_path = os.path.join(output_dir, f"{dataset}-final.png")
        fig.savefig(combined_path, dpi=300, bbox_inches="tight")
        combined_path = os.path.join(output_dir, f"{dataset}-final.pdf")
        fig.savefig(combined_path, dpi=300, bbox_inches="tight")
        plt.close(fig)

        # 保存单独的子图为 PDF
        fig1.savefig(os.path.join(output_dir, f"{dataset}-fig1.pdf"), dpi=300, bbox_inches="tight")
        fig2.savefig(os.path.join(output_dir, f"{dataset}-fig2.pdf"), dpi=300, bbox_inches="tight")
        fig3.savefig(os.path.join(output_dir, f"{dataset}-fig3.pdf"), dpi=300, bbox_inches="tight")

        print(f"Combined figure saved to {combined_path}")
        print(f"Individual figures saved as:")
        print(f"  - {dataset}-fig1.pdf")
        print(f"  - {dataset}-fig2.pdf")
        print(f"  - {dataset}-fig3.pdf")


# In[161]:


if use_spefic_color:
    combine_three_plots(fig1, fig2, fig3, dataset=dataset)


# In[ ]:




