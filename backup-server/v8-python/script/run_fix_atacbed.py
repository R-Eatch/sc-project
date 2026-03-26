#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd


# In[ ]:


def fix_bed(file_path,output_path):
        # 定义列名
    columns = [
        "chrom", "chromStart", "chromEnd", "name", "score", 
        "strand", "signalValue", "pValue", "qValue", "peak"
    ]
    
    # 读取文件（假设是以制表符分隔）
    df = pd.read_csv(file_path, sep="\t", header=None, names=columns)
    
    # 确保 score 和 peak 是整数
    df["score"] = df["score"].round().astype(int)
    df["peak"] = df["peak"].round().astype(int)
    df["signalValue"] = df["signalValue"].round().astype(int)
    df["pValue"] = df["pValue"].round().astype(int)
    df["qValue"] = df["qValue"].round().astype(int)
    # 输出修复后的文件
    df.to_csv(output_path, sep="\t", header=False, index=False)
    print(f"修复后的文件已保存为 {output_path}")
    


# In[ ]:


# 读取 narrowPeak 文件
file_path = "J103660_R-Ges23-MG-bATAC_peaks.narrowPeak"
output_path = "J103660_R-Ges23-MG-bATAC_peaks_fixed.narrowPeak"


# In[6]:


fix_bed(file_path = "J103660_R-Ges23-MG-bATAC_peaks.narrowPeak",output_path = "J103660_R-Ges23-MG-bATAC_peaks_fixed.narrowPeak")


# In[6]:


fix_bed(file_path = "J103660_S-GES14-MG-bATAC_peaks.narrowPeak",output_path = "J103660_S-GES14-MG-bATAC_peaks_fixed.narrowPeak")


# In[ ]:




