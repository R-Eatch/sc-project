import scanpy as  sc
import anndata as ad
#ea = ad.read_h5ad('./processed/allsc_sg_v9_scrublet.h5ad')
ea = ad.read_h5ad('./H-MG.h5ad')
if ea.raw is not None and ea.raw.X is not None:
    print('have raw data')
    adata = ad.AnnData(X=ea.raw.X.copy(),
                       var=ea.raw.var.copy(),
                       obs=ea.obs.copy())
else:
    print('no raw data')
sc.settings.figdir = ''
# 查看表达矩阵 X 的形状
print("X shape:", ea.shape)

# 查看细胞元数据（obs）的前5行
print("\n=== Observations (obs) head ===")
print(ea.obs.head())

# 查看基因元数据（var）的前5行
print("\n=== Variables (var) head ===")
print(ea.var.head())

# 查看 AnnData 对象的整体结构
print("\n=== AnnData object summary ===")
print(ea)
sc.settings.figdir = ''
sc.pl.umap(
        ea,
        color=["author_cell_type", "cell_type",'development_stage','cell_state'],
        # increase horizontal space between panels
        wspace=0.5,
        size=3,
    ncols=3,
    save='H-MG_celltype.png',
    color_map='viridis'
    )

import anndata as ad
import numpy as np
import os


num_cells_to_extract = 20000 
input_file = 'H-MG.h5ad'
output_file = 'H-MG_subset_20k.h5ad'
my_seed = 2025 
ad.read_h5ad(input_file)
def subsample_save_adata(
    adata: ad.AnnData,
    n_cells_to_sample: int,
    output_h5ad_path: str,
    random_seed: int = None
) -> None:
    """
    随机抽取指定数量的细胞 (observations) 从一个 AnnData 对象中，
    并将结果子集保存为一个新的 H5AD 文件。

    Args:
        adata: 输入的 AnnData 对象。
        n_cells_to_sample: 需要随机抽取的细胞数量。
        output_h5ad_path: 子集 AnnData 对象将要保存的文件路径
                          (例如: 'path/to/subset.h5ad')。
        random_seed: 可选整数，用于设置 numpy 随机数生成器的种子，
                     以确保可重复的抽样。如果为 None，则每次抽样结果不同。

    Returns:
        None. 该函数将子集 AnnData 对象保存到磁盘。

    Raises:
        ValueError: 如果 n_cells_to_sample 大于输入 AnnData 对象中的总细胞数。
        TypeError: 如果输入的 adata 不是 AnnData 对象，或者参数类型不正确。
    """
    # --- 输入验证 ---
    if not isinstance(adata, ad.AnnData):
        raise TypeError("输入 'adata' 必须是一个 AnnData 对象。")
    if not isinstance(n_cells_to_sample, int) or n_cells_to_sample <= 0:
         raise ValueError("'n_cells_to_sample' 必须是一个正整数。")
    if not isinstance(output_h5ad_path, str):
        raise TypeError("'output_h5ad_path' 必须是一个字符串。")

    total_cells = adata.n_obs
    print(f"输入的 AnnData 对象包含 {total_cells} 个细胞。")

    # --- 检查请求的细胞数是否合理 ---
    if n_cells_to_sample > total_cells:
        raise ValueError(
            f"请求抽取 {n_cells_to_sample} 个细胞，但输入的 AnnData "
            f"对象中只有 {total_cells} 个可用细胞。"
        )

    # --- 设置随机种子以保证可重复性 ---
    if random_seed is not None:
        print(f"为抽样设置随机种子为: {random_seed}")
        np.random.seed(random_seed)

    print(f"正在随机抽取 {n_cells_to_sample} 个细胞 (无放回)...")

    # --- 生成要保留的细胞的随机索引 ---
    # np.arange(total_cells) 创建一个包含 [0, 1, ..., total_cells-1] 的数组
    # np.random.choice 从中无放回地选择 n_cells_to_sample 个索引
    random_indices = np.random.choice(
        np.arange(total_cells),
        size=n_cells_to_sample,
        replace=False  # 确保是无放回抽样
    )

    # --- 创建子集 AnnData 对象 ---
    # 使用随机索引对 AnnData 对象进行切片。
    # adata[indices, :] 选择指定索引的行（细胞）和所有列（基因/特征）。
    # .copy() 确保我们得到一个独立的对象，而不是原始数据的视图。
    subset_adata = adata[random_indices, :].copy()

    print(f"已创建子集 AnnData，包含 {subset_adata.n_obs} 个细胞和 {subset_adata.n_vars} 个变量。")
    # --- 保存子集 AnnData 对象 ---
    print(f"正在将子集保存到: {output_h5ad_path}")
    try:
        # 使用 gzip 压缩以节省空间
        subset_adata.write_h5ad(output_h5ad_path)
        print("子集成功保存。")
    except Exception as e:
        print(f"保存文件时出错: {e}")
        # 根据需要，您可以选择重新抛出异常
        # raise e

# --- 如何使用这个函数 ---
if __name__ == '__main__':

    # --- 3. 调用函数 ---
    try:
        subsample_save_adata(
            adata=ad.read_h5ad(input_file),
            n_cells_to_sample=num_cells_to_extract,
            output_h5ad_path=output_file,
            random_seed=my_seed
        )

        # --- 4. (可选) 验证输出文件 ---
        print(f"\n验证保存的文件: {output_file}")
        if os.path.exists(output_file):
            loaded_subset = ad.read_h5ad(output_file)
            print("加载的子集信息:")
            print(loaded_subset)
            assert loaded_subset.n_obs == num_cells_to_extract, "抽取的细胞数量不匹配!"
            print("验证成功。")
            # # 如果需要，可以取消下面一行的注释来删除生成的测试文件
            # os.remove(output_file)
            # print(f"已清理测试文件: {output_file}")
        else:
            print(f"错误: 输出文件 {output_file} 在保存后未找到。")

    except ValueError as ve:
        print(f"捕获到预期的 ValueError: {ve}")
    except TypeError as te:
        print(f"捕获到预期的 TypeError: {te}")
    except Exception as e:
        print(f"发生意外错误: {e}")

    # --- 示例：触发 ValueError ---
    print("\n--- 尝试抽取过多细胞 (预期会报错) ---")
    try:
        subsample_save_adata(dummy_adata, n_original_cells + 1, "too_many.h5ad")
    except ValueError as ve:
        print(f"成功捕获到预期的 ValueError: {ve}")
