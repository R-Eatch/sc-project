import os
import anndata as ad
import scanpy as sc

def split_h5ad_file(input_file):
    # Read the input h5ad file
    adata = ad.read_h5ad(input_file)
    
    # 1) Create a new/overwrite species column from the first letter of 'sample'
    adata.obs['species'] = adata.obs['sample'].str[0]

    # 2) Define conditions for splitting based on species and gland
    #    (including the two new groups 'S-SG' and 'M-SG')
    conditions = {
        "M-MG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "MG"),
        "R-MG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "MG"),
        "R-AG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "AG"),
        "S-MG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "MG"),
        "S-AG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "AG"),
        "R-CG": (adata.obs['species'] == "R") & (adata.obs['gland'] == "CG"),
        "S-SG": (adata.obs['species'] == "S") & (adata.obs['gland'] == "SG"),
        "M-SG": (adata.obs['species'] == "M") & (adata.obs['gland'] == "SG")
    }

    # Base directory for the output files
    base_dir = "../v8-python"

    # Process each condition and save the subset
    for group_name, condition in conditions.items():
        subset_adata = adata[condition].copy()
        group_path = os.path.join(base_dir, group_name)
        if not os.path.exists(group_path):
            os.makedirs(group_path)
        subset_adata.write(os.path.join(group_path, f"{group_name}.h5ad"))
        print(f"Saved {group_name}.h5ad to {group_path}")

if __name__ == "__main__":
    # Define the input file path directly in the script
    input_file_path = "./processed/v9_all_counts.h5ad"

    # Call the function with the direct input file path
    split_h5ad_file(input_file_path)

