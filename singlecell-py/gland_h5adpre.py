import scanpy as sc
import os
import gc

# 你的原始文件列表
files = ['./data/MG_pp.h5ad', './data/SG_pp.h5ad', './data/EG_pp.h5ad']
target_n_cells = 50000  # 设定抽样目标

for file_path in files:
    print(f"========== 正在预处理: {file_path} ==========")
    
    # 1. 读取数据
    try:
        adata = sc.read_h5ad(file_path)
    except Exception as e:
        print(f"读取失败: {e}")
        continue
        
    print(f"  原始维度: {adata.shape}")

    # 2. 【核心】将 X 重置为 Counts
    # Sceasy 默认会把 .X 里的数据转成 Seurat 的 counts/data
    if 'counts' in adata.layers:
        print("  - 检测到 layers['counts']，正在替换 adata.X ...")
        adata.X = adata.layers['counts'].copy()
    else:
        print("  - 警告: 未检测到 layers['counts']，将直接使用当前 adata.X")

    # 3. 【核心】随机抽样 (Subsampling)
    if adata.n_obs > target_n_cells:
        print(f"  - 正在抽样: {adata.n_obs} -> {target_n_cells}")
        sc.pp.subsample(adata, n_obs=target_n_cells, random_state=42)
    else:
        print("  - 细胞数不足 50k，跳过抽样")

    # 4. 【关键】清理不必要的 slots，防止 sceasy 报错
    # sceasy 转换时，如果 raw 或某些 unstructured 数据格式不对容易报错
    # 为了 CellChat，我们只需要 X, obs, var
    if adata.raw is not None:
        print("  - 移除 .raw (避免转换混淆)")
        del adata.raw
    
    # 清空 uns 中可能导致 R 读取失败的大型对象 (可选，视情况而定)
    # adata.uns = {} 

    # 5. 保存为临时文件
    # 文件名加上 _subset 后缀
    sample_name = os.path.basename(file_path).replace('.h5ad', '')
    out_path = f"./data/{sample_name}_subset_for_sceasy.h5ad"
    
    print(f"  - 保存预处理后的文件: {out_path}")
    adata.write_h5ad(out_path)
    
    # 强制垃圾回收，释放内存
    del adata
    gc.collect()

print("\n所有预处理完成！现在可以使用 sceasy 转换生成的 _subset 文件了。")
