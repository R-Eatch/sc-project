#!/usr/bin/env python
# coding: utf-8

import scanpy as sc
import anndata as ad
import os

# ==============================================================================
# 1. 配置区域
# ==============================================================================

# 输入文件的根目录
INPUT_ROOT_DIR = '../' 

# 输出目录
OUTPUT_DIR = './data'
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

# 【路径模板】
# python 的 f-string 格式，{prefix} 会被替换为数据集名称 (如 M-MG)
# 假设文件路径是: ../M-MG/M-MG_for_subset.h5ad
FILE_PATH_TEMPLATE = os.path.join(INPUT_ROOT_DIR, '{prefix}','0.loom', '{prefix}_for_subset.h5ad')

# 【分组与重命名规则】
# Key = 新的腺体名称 (Target Group)
# Value = 包含的原始数据集前缀列表
DATASET_GROUPS = {
    # 1. MG 组 (保持 MG)
    'MG': ['M-MG', 'R-MG', 'S-MG'],
    
    # 2. EG 组 (原来的 SG 改名为 EG)
    'EG': ['M-SG', 'S-SG'], 
    
    # 3. SG 组 (原来的 AG 和 CG 合并为 SG)
    'SG': ['R-AG', 'S-AG', 'R-CG'] 
}

# Metadata 中存储腺体信息的列名 (会自动修改此列)
GLAND_COL = 'gland'

# ==============================================================================
# 2. 执行逻辑
# ==============================================================================

def process_and_merge_group(target_group, prefixes):
    print(f"\n{'='*50}")
    print(f"Processing Target Group: {target_group}")
    print(f"{'='*50}")
    
    adatas = []
    
    for prefix in prefixes:
        # 构建具体文件路径
        fpath = FILE_PATH_TEMPLATE.format(prefix=prefix)
        
        if not os.path.exists(fpath):
            print(f"  [MISSING] File not found: {fpath}")
            continue
            
        print(f"  Loading: {prefix} ...")
        try:
            # 读取 h5ad
            adata = sc.read_h5ad(fpath)
            
            # --- 关键步骤：修改腺体标签 ---
            # 无论原来叫什么，现在统一改为 target_group (MG, EG, 或 SG)
            adata.obs[GLAND_COL] = target_group
            
            # 记录原始数据集来源，方便后续分析
            adata.obs['batch_source'] = prefix
            
            # 确保 index 唯一 (加上前缀)
            # 例如: M-MG-AACTGC...
            #adata.obs_names = f"{prefix}-" + adata.obs_names
            
            print(f"    -> Loaded shape: {adata.shape}, Gland set to: '{target_group}'")
            adatas.append(adata)
            
        except Exception as e:
            print(f"    [ERROR] Failed to load {fpath}: {e}")
            
    # --- 合并 ---
    if not adatas:
        print(f"  [WARN] No valid data for {target_group}. Skipping.")
        return

    print(f"  Merging {len(adatas)} datasets...")
    
    # 使用 outer join 合并
    # 即使不同物种/批次的基因列表略有不同，outer join 也会保留并集
    # fill_value=0 确保缺失的表达量填为 0
    adata_merged = ad.concat(
        adatas, 
        join='outer', 
        merge='unique', # 保留 obs 中不重复的列
        fill_value=0
    )
    
    # 再次去重 var_names (防止不同文件中有重复定义的基因名)
    #adata_merged.var_names_make_unique()
    
    print(f"  [SUCCESS] Merged Shape: {adata_merged.shape}")
    print(f"  Example obs check: {adata_merged.obs[[GLAND_COL, 'batch_source']].head(2)}")
    
    # --- 保存 ---
    out_name = f"Final_{target_group}_merged.h5ad"
    save_path = os.path.join(OUTPUT_DIR, out_name)
    
    print(f"  Saving to {save_path} ...")
    adata_merged.write(save_path)
    print("  Done.")

# ==============================================================================
# 3. 主程序
# ==============================================================================

if __name__ == "__main__":
    for group_name, prefix_list in DATASET_GROUPS.items():
        process_and_merge_group(group_name, prefix_list)
        
    print("\nAll merge tasks completed.")
