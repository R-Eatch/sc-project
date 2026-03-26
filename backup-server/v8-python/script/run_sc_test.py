#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scanpy as  sc
import anndata as ad
import loompy
import scvelo as scv
import scanpy as sc
import os
import pickle
import pandas as pd
print('import successful')


# In[2]:


scv.logging.print_version()

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization


# In[2]:


###global variable###
PCs = 10 
res = 0.2
dataset = 'R-MG'
Featuregenes = ['Esr1','Epcam','Top2a']
###########################


# In[46]:


adata =sc.read_h5ad("D:/111/S-AG_cleaned.h5ad")


# In[47]:


adata.X=adata.raw.X


# In[48]:


adata.obs['stage']


# In[49]:


gene_expression = adata[:, 'Prlr'].X
gene_expression


# In[50]:


gene_expression_dense = gene_expression.toarray().flatten()
gene_expression_dense


# In[51]:


# 获取stage信息
stages = adata.obs['stage']

# 将基因表达数据和stage信息合并为DataFrame
expression_df = pd.DataFrame({
    'expression': gene_expression_dense,
    'stage': stages
})

# 按stage分组，计算每个stage的平均表达
mean_expression_by_stage = expression_df.groupby('stage')['expression'].mean()

# 显示结果
print(mean_expression_by_stage)


# In[52]:


import seaborn as sns
import matplotlib.pyplot as plt

# 创建小提琴图
plt.figure(figsize=(8, 6))
sns.violinplot(x='stage', y='expression', data=expression_df)
plt.title('GENE_X Expression Across Stages')
plt.show()


# In[16]:


sc.pl.umap(
        adata,
        color=['celltype',*Featuregenes,'stage','Acta2','sample'],
        wspace=0.5,
        size=3,
    color_map='viridis'
    )


# In[22]:


import matplotlib.pyplot as plt
import seaborn as sns

# 绘制未经过log1p处理的数据分布
plt.figure(figsize=(12, 6))

# 选取前几个基因（特征）进行查看
sns.histplot(adata.X[:, 0].toarray(), kde=True, bins=50)
plt.title("Distribution of the first feature (before log1p?)")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.show()


# In[4]:


adata_raw = sc.read_10x_h5("D:/111/filtered_feature_bc_matrix.h5")


# In[5]:


loom_adata = ad.read_loom("D:/111/S-3MTh-MG.loom")


# In[6]:


df=pd.read_csv("gene_map.csv")


# In[7]:


df


# In[8]:


adata_raw.var


# In[9]:


loom_adata.var


# In[10]:


reference_to_gene=dict(zip(df['reference_id'],df['gene_name']))


# In[11]:


for key,value in reference_to_gene.items():
    if value == '-':
        reference_to_gene[key]=key


# In[12]:


adata_raw.var['Gene']=adata_raw.var['gene_ids'].map(reference_to_gene)


# In[13]:


adata_raw.var.set_index('Gene',inplace=True)


# In[14]:


loom_adata.obs_names=loom_adata.obs_names.str.replace('x','',regex=False)
loom_adata.obs_names


# In[15]:


adata_raw.obs_names.name= 'CellID'
adata_raw.obs_names=adata_raw.obs_names.str.replace('-1','',regex=False)
adata_raw.obs_names='S-3MTh-MG:'+adata_raw.obs_names
adata_raw.obs_names


# In[16]:


adata_raw.obs['uid']=adata_raw.obs.index


# In[17]:


common_cells = adata_raw.obs_names.intersection(loom_adata.obs_names)
common_genes = adata_raw.var_names.intersection(loom_adata.var_names)
adata_sub = adata_raw[adata_raw.obs_names.isin(common_cells),adata_raw.var_names.isin(common_genes)].copy()
ldata_sub = loom_adata[loom_adata.obs_names.isin(common_cells),loom_adata.var_names.isin(common_genes)].copy()
print(f'cells:{len(common_cells)},genes:{len(common_genes)} {"#"*80}',flush=True)
adata = scv.utils.merge(adata_sub, ldata_sub)
adata.write('S-3MTh-MG_for_subset.h5ad')


# In[18]:


####Quality Control(temporary use default)####


# In[19]:


adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")


# In[20]:


sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)


# In[21]:


sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)


# In[22]:


sc.pp.scrublet(adata)


# In[23]:


adata.raw = adata
adata.layers["counts"] = adata.X.copy()


# In[24]:


# Normalizing to median total counts
sc.pp.normalize_total(adata,target_sum= 1e4)
# Logarithmize the data
sc.pp.log1p(adata)
print('Finish normalized')


# In[25]:


sc.pp.highly_variable_genes(adata, n_top_genes=3000)
#ea = ea[:, ea.var.highly_variable]
print('Finish Varible genes')


# In[26]:


sc.pp.scale(adata,max_value=10)
sc.tl.pca(adata,n_comps=PCs)
sc.pl.pca_variance_ratio(adata,n_pcs=PCs,log = True)


# In[27]:


sc.pp.neighbors(adata,n_pcs= PCs)
#sc.tl.louvain(ea,resolution = res)
sc.tl.leiden(adata,resolution = res)
print('Finish clustering')
sc.tl.umap(adata)
print('Finish UMAP')


# In[47]:


sc.pl.umap(
        adata,
        color=['leiden',*Featuregenes],
        wspace=0.5,
        size=3,
    color_map='viridis'
    )


# In[29]:


adata.X=adata.raw.X


# In[30]:


scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000,enforce=True)
sc.pp.neighbors(adata, n_pcs=10, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)


# In[31]:


scv.tl.velocity(adata)


# In[32]:


scv.tl.velocity_graph(adata)


# In[33]:


scv.pl.proportions(adata)


# In[34]:


scv.pl.velocity_embedding_stream(adata, basis='umap',color='leiden')


# In[35]:


adata.obs


# In[36]:


adata.var


# In[37]:


adata


# In[ ]:




