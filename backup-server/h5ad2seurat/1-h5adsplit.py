import anndata as ad
import scipy.sparse as sp
import os

# 分割 h5ad 文件
def split_h5ad_file(input_file):
    adata = ad.read_h5ad(input_file)
    conditions = {
        "M-MG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "MG"),
        "R-MG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "MG"),
        "R-AG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "AG"),
        "S-MG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "MG"),
        "S-AG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "AG"),
        "R-CG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "CG")
    }

    output_dir = "../newsubset"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for group_name, condition in conditions.items():
        subset_adata = adata[condition].copy()
        group_path = os.path.join(output_dir, group_name)
        if not os.path.exists(group_path):
            os.makedirs(group_path)
        subset_adata.write(os.path.join(group_path, f"{group_name}.h5ad"))
        print(f"Saved {group_name}.h5ad to {group_path}")

# 提取原始表达矩阵和元数据
def extract_raw_data(input_dir):
    output_dir = "../rawsubset"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    conditions = ["M-MG", "R-MG", "R-AG", "S-MG", "S-AG", "R-CG"]
    for condition in conditions:
        input_file_path = os.path.join(input_dir, condition, f"{condition}.h5ad")
        output_file_path = os.path.join(output_dir, f"{condition}_raw.h5ad")
        
        adata = ad.read_h5ad(input_file_path)
        if raw_data is not None:
            raw_matrix = raw_data.X
            raw_obs = raw_data.obs
            raw_var = raw_data.var
        else:
            raw_matrix = adata.X
            raw_obs = adata.obs
            raw_var = adata.var
        
        new_adata = ad.AnnData(X=raw_matrix, obs=raw_obs, var=raw_var)
        new_adata.write_h5ad(output_file_path)
        print(f"Saved raw data of {condition} to {output_file_path}")

if __name__ == "__main__":
    input_file = "allcell_3species_v5_noSG.h5ad"
    print("Splitting h5ad file...")
    split_h5ad_file(input_file)
    print("Splitting completed.")
    
    print("Extracting raw data...")
    extract_raw_data("../newsubset")
    print("Extraction completed.")

