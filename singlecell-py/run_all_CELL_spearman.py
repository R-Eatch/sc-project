#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import scipy
import scanpy as  sc
import anndata as ad
import warnings
import seaborn as sns
warnings.filterwarnings("ignore", category=UserWarning, message="ImplicitModificationWarning")
print('import successful')


# In[2]:


####global variable###


# In[3]:


datasets = ['M-MG','R-MG','S-MG','R-AG','S-AG','R-CG']
use_hvg=True
n_top_genes=2000
adata_list=[]
################


# In[ ]:


adata_combined=sc.read_h5ad("../../h5ad2seurat/allspecies_v8_nosg_raw_re2.h5ad")


# In[6]:


adata_combined.obs['celltype_dataset']= (adata_combined.obs['celltype'].astype(str)
                                         + "_" + adata_combined.obs['species'].astype(str)
                                         + "_" + adata_combined.obs['gland'].astype(str))


# In[ ]:


if use_hvg:
    sc.pp.highly_variable_genes(adata_combined, n_top_genes=n_top_genes)
    adata_combined = adata_combined[:, adata_combined.var["highly_variable"]].copy()


# In[8]:


celltype_datasets=adata_combined.obs['celltype_dataset'].unique()
celltype_datasets


# In[19]:


celltype_subsets = {
    celltype: adata_combined[adata_combined.obs['celltype_dataset'] == celltype]
    for celltype in celltype_datasets
}

# 存储结果
correlation_results = []

for celltype1 in celltype_datasets:
    for celltype2 in celltype_datasets: 
        subset_gland1 = celltype_subsets[celltype1]
        subset_gland2 = celltype_subsets[celltype2]

        if subset_gland1.shape[0] <= 1 or subset_gland2.shape[0] <= 1:
            print(f"{celltype1} or {celltype2} has too few samples, skipping correlation.")
            continue

        mean_expression_gland1 = np.asarray(subset_gland1.X.mean(axis=0)).flatten()
        mean_expression_gland2 = np.asarray(subset_gland2.X.mean(axis=0)).flatten()

        subset_gland1.X = np.nan_to_num(subset_gland1.X, nan=mean_expression_gland1)
        subset_gland2.X = np.nan_to_num(subset_gland2.X, nan=mean_expression_gland2)

        mean_expression_gland1 = np.asarray(subset_gland1.X.mean(axis=0)).flatten()
        mean_expression_gland2 = np.asarray(subset_gland2.X.mean(axis=0)).flatten()

        corr, pval = spearmanr(mean_expression_gland1, mean_expression_gland2)

        # 保存结果
        correlation_results.append({
            'celltype1': celltype1,
            'celltype2': celltype2,
            'correlation': corr,
            'p-value': pval
        })
        print(f'{celltype1} vs {celltype2} correlation: {corr}, p-value: {pval}')


correlation_df = pd.DataFrame(correlation_results).dropna()
correlation_df = correlation_df[correlation_df['p-value'] < 0.05]
print(correlation_df)


# In[18]:


enumerate(celltype_datasets)


# In[20]:


correlation_matrix = correlation_df.pivot(index='celltype1', columns='celltype2', values='correlation')
correlation_matrix.to_csv('ALL-correlation_matrix.csv')


# In[21]:


plt.figure(figsize=(20,18))  # 调整图像大小
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', vmin=-1, vmax=1, 
            cbar_kws={'label': 'Spearman Correlation'}, linewidths=0.5)
plt.title('Spearman Correlation between Celltype1 and Celltype2', fontsize=16)
plt.xlabel('Celltype2', fontsize=12)
plt.ylabel('Celltype1', fontsize=12)

# 调整坐标轴标签的旋转角度和大小
plt.xticks(rotation=90, fontsize=8)
plt.yticks(rotation=0, fontsize=8)

# 保存图像
plt.savefig('ALL-heatmap_correlation.png', dpi=300, bbox_inches='tight')  # 保存图像
plt.show()  # 显示图像


# In[22]:


correlation_matrix1= correlation_matrix.fillna(0)
# 使用 Seaborn 绘制带聚类的热图
plt.figure(figsize=(20,18))  # 调整图像大小
sns.clustermap(correlation_matrix1, annot=True, cmap='coolwarm', vmin=-1, vmax=1, 
               figsize=(16, 14), cbar_kws={'label': 'Spearman Correlation'}, 
               linewidths=0.5, row_cluster=True, col_cluster=True)

# 设置标题和保存图像
plt.suptitle('Spearman Correlation between Celltype1 and Celltype2 (Clustermap)', fontsize=16)
plt.subplots_adjust(top=0.93)  # 调整标题位置，避免与热图重叠

plt.savefig('ALL-heatmap_correlation-cluster.png', dpi=300, bbox_inches='tight')  # 保存图像
plt.show()  # 显示图像


# In[ ]:




