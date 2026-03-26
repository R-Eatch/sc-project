#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd

df = pd.read_csv("D:/111/rabbit.gtf", sep='\t', comment='#', header=None)

df1=df[df[2]=='gene']

genes_df=df[df[2]=='gene']

genes_df['gene_id']=df1[8].str.extract('gene_id "([^"]+)"')
genes_df['gene_name']=df1[8].str.extract('gene_name "([^"]+)"')

genes_df[[0,3,4,'gene_name',6]]

bed_df=genes_df[[0,3,4,'gene_name',6]]
bed_df[3]=bed_df[3]-1
bed_df

alias_txt=genes_df[['gene_id','gene_name']]

alias_txt['alias'] = alias_txt['gene_name'] + '&' + alias_txt['gene_name']

bed_df.to_csv('genes_regions.bed',sep='\t',header=None,index=False)
alias_txt.to_csv('genes_alias.txt',sep='\t',header=None,index=False)


