#!/usr/bin/env python
# coding: utf-8

import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import plotly.express as px
import matplotlib.pyplot as plt
import matplotlib
import scvelo as scv
import os

print('import successful')

# --- Scvelo Settings ---
scv.logging.print_version()
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.settings.figdir = './figure_velo'  # 设置统一的保存文件夹
os.makedirs(scv.settings.figdir, exist_ok=True)

# --- Global Variables ---
# 对应 gland_pp.py 的列表
glandlist = ['MG', 'SG', 'EG'] 
use_dynamic_method = False
run_score = True
PCs = 50  # 建议与 gland_pp.py 保持一致
n_neighbors = 15 # 建议与 gland_pp.py 保持一致
n_top_genes = 3000

# 如果 pp 阶段已经做过 clean，这里通常不需要再 drop，除非有特定需求
drop_cells = False 
droplist = ['Lum-Cd74', 'Lum-Adipoq']

# --- Functions ---

def run_plotly(adata, dataset):
    # 检查是否存在 UMAP，防止报错
    if 'X_umap' not in adata.obsm:
        print(f"Skipping plotly: X_umap not found in {dataset}")
        return

    umap_data = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'])
    umap_data['cell_ids'] = adata.obs.index 
    # 确保 metadata 存在
    if 'stage' in adata.obs:
        umap_data['stage'] = adata.obs['stage']
        hover_data = ['cell_ids', 'stage']
    else:
        hover_data = ['cell_ids']

    fig = px.scatter(umap_data, x='UMAP1', y='UMAP2', hover_data=hover_data, width=800, height=450)
    # fig.show() # 在服务器/脚本模式下通常不 show，直接保存
    fig.write_html(f'./figure_velo/{dataset}-interactive_umap_plot.html')
    print(f"Saved interactive plot for {dataset}")


def process_velocity(dataset):
    print(f'--- Processing {dataset} ---')
    
    # 1. 读取 gland_pp.py 的输出文件
    file_path = f'./data/{dataset}_pp.h5ad'
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}, skipping...")
        return
        
    adata = sc.read_h5ad(file_path)
    
    # 2. (可选) 二次过滤细胞
    if drop_cells and 'newcelltype' in adata.obs:
        adata = adata[~adata.obs['newcelltype'].isin(droplist)].copy()
        print(f"Dropped cells: {droplist}")
    
    # 3. 绘制交互式 Plotly (基于已有 UMAP)
    run_plotly(adata, dataset)
    
    # 4. scVelo 预处理
    # 注意：filter_and_normalize 主要针对 spliced/unspliced layers 进行标准化
    # enforce=True 强制保留之前选定的 highly_variable_genes (如果需要严格对齐)
    scv.pp.filter_and_normalize(adata, n_top_genes=n_top_genes)
    
    # 5. 关键步骤：计算 Neighbors 和 Moments
    # 为了保持与 gland_pp.py 的聚类一致性，必须使用相同的 PCA 表示
    if 'X_pca_harmony' in adata.obsm:
        print("Using Harmony embeddings for velocity neighbors...")
        use_rep = 'X_pca_harmony'
    elif 'X_pca' in adata.obsm:
        print("Using PCA embeddings for velocity neighbors...")
        use_rep = 'X_pca'
    else:
        print("Warning: No PCA/Harmony found, recalculating PCA (might differ from pp)...")
        sc.tl.pca(adata)
        use_rep = 'X_pca'

    # 计算 neighbors (用于 moments 平滑)
    sc.pp.neighbors(adata, n_pcs=PCs, n_neighbors=n_neighbors, use_rep=use_rep, random_state=42)
    
    # 计算 moments (velocity 的基础)
    scv.pp.moments(adata, n_pcs=PCs, n_neighbors=n_neighbors)
    
    # 6. 计算 Velocity
    if use_dynamic_method:    
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
        print("Mode: dynamical")
    else:
        scv.tl.velocity(adata, mode='stochastic')
        print("Mode: stochastic (default)")
        
    scv.tl.velocity_graph(adata)
    
    # 7. 绘图与保存
    # 确定绘图用的分类标签
    color_key = 'celltype' # 默认为 celltype，如果你的 pp 脚本里是 newcelltype 请修改
    #if 'newcelltype' in adata.obs:
    #    color_key = 'newcelltype'
    #elif 'leiden' in adata.obs:
    #    color_key = 'leiden'
        
    print(f"Plotting with color: {color_key}")

    # Stream Plot (PNG)
    scv.pl.velocity_embedding_stream(
        adata, basis='umap', color=color_key, 
        save=f'{dataset}_velo.png', title=f'{dataset} Velocity'
    )
    
    # Stream Plot (PDF with Rasterization)
    ax = scv.pl.velocity_embedding_stream(
        adata, basis='umap', color=color_key, show=False, title=f'{dataset} Velocity'
    )
    if ax is not None: # 防止绘图失败
        fig = ax.get_figure()
        for artist in fig.findobj(match=matplotlib.artist.Artist):
            if hasattr(artist, "set_rasterized"):
                artist.set_rasterized(True)
        fig.savefig(f'{scv.settings.figdir}/{dataset}_velo.pdf', 
                    dpi=300, bbox_inches="tight", format="pdf")
        plt.close(fig)

    # Proportions Plot
    scv.pl.proportions(adata, groupby=color_key, save=f'{dataset}_proportions.png')
    if 'leiden' in adata.obs:
        scv.pl.proportions(adata, groupby='leiden', save=f'{dataset}_proportions-leiden.png')

    # 8. (可选) Cell Cycle Score
    if run_score:
        # 只有当基因名匹配时才有效
        try:
            scv.tl.score_genes_cell_cycle(adata)
            scv.pl.scatter(adata, 
                        color_gradients=['S_score', 'G2M_score'], 
                        palette=['green', 'orange'], 
                        smooth=True, perc=[5, 90],
                        save=f'{dataset}_cycle.png'
                        )
        except Exception as e:
            print(f"Skipping cell cycle scoring: {e}")

    # 9. Pseudotime
    try:
        scv.tl.velocity_pseudotime(adata)
        scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save=f'{dataset}_velo_pseu.png')
    except Exception as e:
        print(f"Skipping pseudotime: {e}")
    
    # 10. 保存最终结果
    adata.write(f'./data/{dataset}_velofinished.h5ad')
    print(f"{dataset} finished and saved.")


# --- Main Execution ---
def main():
    for dataset in glandlist:
        process_velocity(dataset)

if __name__ == "__main__":
    main()
