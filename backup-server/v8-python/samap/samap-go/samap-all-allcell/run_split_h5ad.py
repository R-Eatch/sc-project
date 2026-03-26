import os
import anndata as ad
import scanpy as sc

def split_h5ad_file(input_file):
    # 读取输入的 h5ad 文件
    adata = ad.read_h5ad(input_file)
    adata.obs['species'] = adata.obs['sample'].str[0]
    # 根据 species 和 gland 定义分割条件
    conditions = {
        "M-MG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "MG"),
        "R-MG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "MG"),
        "R-AG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "AG"),
        "S-MG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "MG"),
        "S-AG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "AG"),
        "R-CG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "CG"),
	"M-SG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "SG"),
	"S-SG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "SG"),

    }

    # 输出文件的基本目录
    base_dir = "./split"
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)

    # 对每个条件进行处理并保存子集到指定目录（不创建额外子文件夹）
    for group_name, condition in conditions.items():
        subset_adata = adata[condition].copy()
        output_file = os.path.join(base_dir, f"{group_name}_ori.h5ad")
        subset_adata.write(output_file)
        print(f"Saved {group_name}.h5ad to {base_dir}")

if __name__ == "__main__":
    # 直接在脚本中定义输入文件路径
    input_file_path = "/data01/sunxuebo/project/scrnaseq/h5ad2seurat/processed/v9_all_counts.h5ad"
    
    # 调用函数，传入输入文件路径
    split_h5ad_file(input_file_path)

