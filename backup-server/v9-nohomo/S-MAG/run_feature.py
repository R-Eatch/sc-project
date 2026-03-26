import scanpy as sc
import matplotlib.pyplot as plt

def plot_feature_genes(
    adata_path,
    feature_genes,
    embedding="X_umap",
    image_format="png",
    dpi=300
):
    """
    从指定的 h5ad 文件读取 AnnData 对象，遍历特征基因列表，对每个存在的基因作图，并输出到对应的图片文件。

    :param adata_path:      你的 AnnData (.h5ad) 文件路径
    :param feature_genes:   需要绘制的基因列表
    :param embedding:       可视化所使用的低维空间（默认为 UMAP）
    :param image_format:    输出图像格式，默认为 png
    :param dpi:             图像分辨率，默认为 300
    """

    # 读取 adata 对象
    adata = sc.read_h5ad(adata_path)
    
    # 根据需要，也可以在这里设置一些 Scanpy 全局参数
    # sc.set_figure_params(dpi_save=300, format='png', ...)
    # 也可以先检查是否存在该 embedding，如 X_umap
    if embedding not in adata.obsm.keys():
        raise ValueError(f"指定的 embedding: {embedding} 不存在于 adata.obsm 中。当前可用: {list(adata.obsm.keys())}")

    for gene in feature_genes:
        # 检查基因名是否在 var_names 中
        if gene in adata.var_names:
            # 绘图
            sc.pl.umap(
                adata,
                color=gene,
                show=False,       # 不弹出窗口，直接保存
                title=gene,       # 标题可设为基因名
                use_raw=False     # 如果需要用 raw 数据，请设置 True，视具体情况而定
            )
            # 保存
            plt.savefig(
                f"{gene}.{image_format}",
                dpi=dpi,
                bbox_inches='tight'
            )
            plt.close()  # 及时关闭图形，避免后续重叠
        else:
            print(f"基因 {gene} 未在 adata.var_names 中找到，跳过。")

    print("绘图完成。")


if __name__ == "__main__":
    # 需要绘制的基因列表
    my_feature_genes = [
        "FGF10","BMP4","WNT6","LEF1","TOP2A","KRT18-1","ESR1","KRT23","KRT14","LALBA","FAAH","PPL","SBSN"
    ]
    
    # 你的 AnnData 文件路径（示例）
    adata_file = "S-MAG_cleaned.h5ad"

    # 调用函数
    plot_feature_genes(adata_file, my_feature_genes)
