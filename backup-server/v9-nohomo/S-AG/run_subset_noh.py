#!/usr/bin/env python
# coding: utf-8

import os
# 固定随机数种子，保证结果可重复（尤其是 UMAP/Leiden 聚类等）
os.environ['PYTHONHASHSEED'] = '0'

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
import seaborn as sns

# 如果需要用到 Scrublet 做双胞检测，需要额外的扩展包：
import scanpy.external as sce


print("Import successful")

###############################
### 全局变量及流程配置示例  ###
###############################

# 全局随机种子
random_seed = 2023
random.seed(random_seed)
np.random.seed(random_seed)

# 输入和输出路径示例
dataset='S-AG'
input_h5ad_path = "/data01/sunxuebo/project/scrnaseq/h5ad2seurat/processed/allsc_sg_v9_scrublet.h5ad"  # <-- 修改成自己的h5ad路径
output_h5ad_path = f"{dataset}_nohomo.h5ad"

# 是否进行 Harmony
run_harmony = True
# 是否进行 scVI
use_scvi = False
# 主成分数量
PCs = 20
# 聚类分辨率
res = 0.3
# 高可变基因数量
top_genes = 2000

# 特征基因（以乳腺分泌细胞相关基因为例，根据需要自行调整）
Featuregenes = [
    "KRT8",
    "KRT18-1",
    "KRT14",
    "ESR1",
    "PGR",
    #"CSN1S1",
    #"CSN2",
    #"CSN3",
    "LALBA",
    "KRT23",
    "FAAH",
    "PLIN2",
    "ACACA",
    "ELF5",
    "EPCAM",
    "PIGR",
    "KRT1",
    "SBSN"
]

###############################
###    函数定义区 (示例)    ###
###############################

def run_qc(adata, doublet_rate=0.06):
    """
    对原始数据进行质量控制并进行双胞检测/去除示例
    其中双胞检测步骤可根据需要自行增减或修改参数
    """
    # 1) 计算线粒体基因百分比
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    
    # 2) 筛掉细胞数或基因数极少的情况（示例阈值，可根据需要修改）
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    
    # 3) 筛掉线粒体含量过高的细胞
    # 假设线粒体含量大于 20% 认为是低质量细胞，可以根据需要调整
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    
    # 4) 双胞检测 (Scrublet)
    print("Running scrublet for doublet detection...")
    sce.pp.scrublet(
        adata,
        expected_doublet_rate=doublet_rate,
        random_state=random_seed
    )
    # 这里会在 adata.obs 中生成 'doublet_score' 和 'predicted_doublet' 等信息
    # 也可以 plot distribution 看看阈值
    # sce.pl.scrublet_score_distribution(adata)
    
    # 5) 去除双胞
    adata = adata[adata.obs["predicted_doublet"] == False].copy()
    
    return adata

def run_preprocess(adata, top_genes=2000):
    """
    进行标准化、对数转换并将结果保存到 layers，随后挑选高可变基因
    """
    # 在对数转换前保存原始counts
    adata.layers["counts"] = adata.X.copy()
    
    # 标准化
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.raw=adata
    # 保存标准化结果
    adata.layers["normalized"] = adata.X.copy()
    print("Finish normalized and store in adata.layers['normalized']")
    
    # 高可变基因
    sc.pp.highly_variable_genes(adata, n_top_genes=top_genes)
    print("Finish selecting HVGs, n_top_genes = ", top_genes)

    return adata

def run_dimension_reduction(adata, use_scvi=False, run_harmony=False, PCs=20, res=0.3):
    """
    包含 PCA、邻近图构建、聚类、UMAP 等流程
    use_scvi / run_harmony 按需修改
    """
    # 对高可变基因进行 scale
    sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, use_highly_variable=True, random_state=random_seed)
    sc.pl.pca_variance_ratio(adata, log=False, show=False,save=f"{dataset}_pcas.png")
    print("Finish PCA")
    
    # 是否使用 scVI，若需要先 setup_anndata
    if use_scvi:
        import scvi
        scvi.settings.seed = random_seed
        print("Running scvi...")
        scvi.model.SCVI.setup_anndata(adata)
        model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
        model.train()
        adata.obsm["X_scVI"] = model.get_latent_representation()
        sc.pp.neighbors(adata, use_rep="X_scVI")
        # 可以额外计算 UMAP
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=res, random_state=random_seed)
        print("Finish scVI + Leiden + UMAP")
        return adata
    else:
        print("Skipping scVI")
    
    # 是否使用 Harmony
    if run_harmony:
        import harmonypy as hm
        print("Running Harmony...")
        pca_result = adata.obsm["X_pca"]
        # 以 adata.obs['sample'] 为例，如果需要多批次矫正，可传入列表
        ho = hm.run_harmony(pca_result, adata.obs, vars_use=["sample"], random_state=random_seed)
        adata.obsm["X_pca_harmony"] = ho.Z_corr.T
        sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=15, n_pcs=PCs, random_state=random_seed)
        print("Finish Harmony + neighbors")
    else:
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=PCs, random_state=random_seed)
        print("Finish neighbors")
    
    # Leiden 聚类
    sc.tl.leiden(adata, resolution=res, random_state=random_seed)
    print("Finish Leiden Clustering")
    
    # UMAP
    sc.tl.umap(adata, random_state=random_seed)
    print("Finish UMAP")
    
    return adata

def plotting(adata, feature_genes=None, save_prefix=""):
    """
    常规可视化示例
    """
    if feature_genes is None:
        feature_genes = []
    
    # 绘制 UMAP，查看聚类结果
    sc.pl.umap(
        adata, 
        color=["sample", "leiden","stage"], 
        wspace=0.5, 
        size=3,
        ncols=2,
        legend_loc="right margin",
        show=False,
        save=f"{save_prefix}_cluster.png"
    )
    
    # 绘制指定特征基因
    sc.pl.umap(
        adata,
        color=feature_genes,
        wspace=0.5,
        size=3,
        ncols=3,
        color_map="plasma",
        show=False,
        save=f"{save_prefix}_features.png"
    )
    print("Finish plotting UMAP & Feature genes")
def export_deg_result(adata):
    """
    Extract DE results from adata.uns['rank_genes_groups'] into a DataFrame.
    Returns only pval < 0.01 subset.
    """
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


def do_deg(adata, sample_prefix="sample"):
    """
    Perform differential expression analysis using rank_genes_groups and save results.
    """
    adata.X=adata.layers['normalized']
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')
    sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby='leiden',
        n_genes=5,
        save=f'{sample_prefix}_dotplot_leiden.png',
        min_logfoldchange=0.25
    )
    df1 = export_deg_result(adata)
    df1.to_csv(f'{sample_prefix}_ranked_genes_leiden.csv', index=False)
    print("DE analysis completed and results saved.")
###############################
###    主流程示例 (main)    ###
###############################

def main():
    # 1. 读取 h5ad
    print(f"Reading AnnData from: {input_h5ad_path}")
    adata = sc.read_h5ad(input_h5ad_path)
    print(adata)
    
    # 2. 根据 adata.obs['sample']，筛选出以“S-AG”开头的样本
    #   假设 sample 格式类似 “S-AG_1”, “S-AG_2” 等
    #   也可能是 “S-AG_01” 等，正则也可。这里用简单 startswith
    mask = [str(s).startswith(f"{dataset}") for s in adata.obs["sample"]]
    adata = adata[mask, :].copy()
    print(f"Subset data by samples starting with {dataset}. Shape now: {adata.shape}")
    adata.write(f"{dataset}_test.h5ad")
    # 3. 进行单细胞数据质控 (包含双胞检测等)
   # adata = run_qc(adata, doublet_rate=0.06)
    print(f"After QC: {adata.shape}")
    adata.obs['sample'] = adata.obs['sample'].replace({'S-AG-1MTH-1': 'S-AG-3MTH-1','S-AG-1MTH-2': 'S-AG-3MTH-2'})
    adata1 = sc.read_h5ad(f"/data01/sunxuebo/project/scrnaseq/v8-python/{dataset}/1.subset/{dataset}_cleaned.h5ad")
    sample_stage_dict = dict(zip(adata1.obs['sample'], adata1.obs['stage']))
    adata.obs['stage'] = adata.obs['sample'].map(sample_stage_dict)
    # 4. 标准化、对数转换，并保存矩阵到 layers；同时挑选高可变基因
    adata = run_preprocess(adata, top_genes=top_genes)
    
    # 5. PCA、邻近图、聚类、UMAP 等
    adata = run_dimension_reduction(
        adata,
        use_scvi=use_scvi,
        run_harmony=run_harmony,
        PCs=PCs,
        res=res
    )
    
    # 6. 绘图：查看聚类结果 + 乳腺分泌细胞特征基因
    plotting(adata, feature_genes=Featuregenes, save_prefix=f"{dataset}_filtered")
    
    # 7. 保存处理完成的 AnnData
    print(f"Writing AnnData to: {output_h5ad_path}")
    adata.write(output_h5ad_path)
    do_deg(adata, sample_prefix=f"{dataset}")
    print("All done!")

if __name__ == "__main__":
    main()

