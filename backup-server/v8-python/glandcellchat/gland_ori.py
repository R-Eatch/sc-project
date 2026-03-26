import scanpy as sc
import pandas as pd
import numpy as np
import os

# ==========================================
# 1. 路径配置
# ==========================================
V9_INPUT_PATH = "../../h5ad2seurat/processed/v9_all_log1p.h5ad"

TARGET_FILES = {
    "MG": "./data/Final_MG_merged.h5ad",
    "SG": "./data/Final_SG_merged.h5ad",
    "EG": "./data/Final_EG_merged.h5ad"
}

# 输出路径配置
PLOT_OUTPUT_DIR = "./output_plots-ori"  # UMAP 可视化 PDF 保存目录
DATA_OUTPUT_DIR = "./data"              # 新产生的 h5ad 数据保存目录

# 自动创建输出文件夹（如果不存在）
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)
os.makedirs(DATA_OUTPUT_DIR, exist_ok=True)

# 设置 Scanpy 的图片输出目录
sc.settings.figdir = PLOT_OUTPUT_DIR

# ==========================================
# 2. 读取总数据并预处理
# ==========================================
print(f"正在加载 V9 总数据集: {V9_INPUT_PATH} ...")
adata_v9 = sc.read_h5ad(V9_INPUT_PATH)

if 'sample' not in adata_v9.obs.columns:
    raise ValueError(f"数据集 {V9_INPUT_PATH} 的 obs 中缺少 'sample' 列，无法重构细胞 ID！")

# 打印 v9 的 obsm 信息
obsm_keys = list(adata_v9.obsm.keys())
print(f"\n[信息] V9 数据集 obsm 中包含以下降维/坐标信息: {obsm_keys}")

# 按照指定逻辑修改 v9 数据的 obs_names
print("正在格式化 V9 数据的细胞 ID (obs_names)...")
adata_v9.obs['cellid'] = (
    adata_v9.obs['sample'].astype(str) + '-' + adata_v9.obs_names.astype(str)
).str.replace(r'(-\d+)+$', '', regex=True)

adata_v9.obs_names = adata_v9.obs['cellid']

# ==========================================
# 3. 提取所有的 obsm 矩阵
# ==========================================
# 为了节省内存并方便映射，我们将所有的 obsm 矩阵提取为带有正确 index 的 Pandas DataFrame 字典
v9_obsm_dfs = {}
for key in obsm_keys:
    # 提取矩阵并去重 index，防止极小概率的正则替换导致细胞名重复报错
    df = pd.DataFrame(adata_v9.obsm[key], index=adata_v9.obs_names)
    v9_obsm_dfs[key] = df[~df.index.duplicated(keep='first')]

# 释放 V9 完整对象以节省内存
del adata_v9 

# ==========================================
# 4. 循环映射 obsm 至目标数据
# ==========================================
for gland, file_name in TARGET_FILES.items():
    if not os.path.exists(file_name):
        print(f"\n[警告] 未找到文件 {file_name}，跳过 {gland} 的处理。")
        continue

    print(f"\n" + "="*50)
    print(f"正在处理 {gland} 腺体数据: {file_name}")
    print("="*50)
    
    adata_target = sc.read_h5ad(file_name)
    target_cells = adata_target.obs_names
    
    # 因为 v9 的细胞名在前面的去重中可能发生微小变化，以提取出的任意一个 DataFrame 的 index 为准
    v9_cells = v9_obsm_dfs[obsm_keys[0]].index

    # 统计匹配情况
    matched_cells = target_cells.intersection(v9_cells)
    missing_in_v9 = target_cells.difference(v9_cells)
    
    print(f"统计信息:")
    print(f"  -> {gland} 细胞总数: {len(target_cells)}")
    print(f"  -> 成功匹配到 v9 的细胞数: {len(matched_cells)}")
    print(f"  -> 缺失 v9 坐标的细胞数: {len(missing_in_v9)}")

    matched_mask = target_cells.isin(matched_cells)

    # 循环遍历并映射所有的 obsm 键
    for key, df in v9_obsm_dfs.items():
        n_dim = df.shape[1] # 获取该降维矩阵的维度 (如 PCA 是 50，UMAP 是 2)
        
        # 创建一个全是 NaN 的空矩阵，形状为 (目标细胞数, 维度)
        new_obsm_arr = np.full((len(target_cells), n_dim), np.nan)
        
        # 提取匹配坐标并赋值
        new_obsm_arr[matched_mask] = df.loc[target_cells[matched_mask]].values
        
        # 赋值回目标 adata 的 obsm，如果存在同名 key 则会被覆盖替换
        adata_target.obsm[key] = new_obsm_arr
        print(f"  -> 成功映射: {key} (维度: {n_dim})")

    # ==========================================
    # 5. 可视化与保存
    # ==========================================
    # 仅当存在 UMAP 时才绘图
    if 'X_umap' in adata_target.obsm.keys():
        print(f"正在为 {gland} 生成并保存 UMAP 图像...")
        adata_plot = adata_target[matched_mask].copy()
        sc.pl.umap(
            adata_plot, 
            title=f"{gland} mapped UMAP (n={len(matched_cells)})",
            save=f"_{gland}_mapped_v9.pdf",
            show=False
        )
    
    # 保存最终数据
    # 使用 os.path.basename 防止拼接出像 ./data/./data/... 这样的冗余路径
    save_filename = os.path.basename(file_name).replace(".h5ad", "_with_v9_obsm.h5ad")
    save_path = os.path.join(DATA_OUTPUT_DIR, save_filename)
    
    print(f"正在保存更新后的数据集至: {save_path} ...")
    adata_target.write_h5ad(save_path)

print("\n全部处理流程执行完毕！")
