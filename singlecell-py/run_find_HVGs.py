#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
from itertools import combinations
import scanpy as sc
import numpy as np
from venn import venn
import os
import warnings

# --- Configuration ---
N_TOP_GENES = 5000
RE_FIND_TOP_GENES = True # Always re-find HVGs based on subsetted data per stage

# Define base file paths and sample names
BASE_DATA_DIR = "../../" # Adjust if needed
SAMPLE_INFO = {
    "R-MG": os.path.join(BASE_DATA_DIR, "R-MG/1.subset/R-MG_cleaned.h5ad"),
    "R-AG": os.path.join(BASE_DATA_DIR, "R-AG/1.subset/R-AG_cleaned.h5ad"),
    "S-MG": os.path.join(BASE_DATA_DIR, "S-MG/1.subset/S-MG_cleaned.h5ad"),
    "S-AG": os.path.join(BASE_DATA_DIR, "S-AG/1.subset/S-AG_cleaned.h5ad"),
    "M-MG": os.path.join(BASE_DATA_DIR, "M-MG/1.subset/M-MG_cleaned.h5ad"),
    "R-CG": os.path.join(BASE_DATA_DIR, "R-CG/1.subset/R-CG_cleaned.h5ad")
}
# 3. 需要提取数据的目标列名（交集名称）
TARGET_INTERSECTIONS = [
    'M-MG & R-MG & S-MG',                     # 三种 MG 的交集
    'M-MG & R-AG & R-CG & R-MG & S-AG & S-MG', # 六个样本的总交集 (请确认此列名与你CSV中的完全一致)
    'R-AG & S-AG'                              # 两种 AG 的交集
]
ALL_SAMPLE_NAMES = list(SAMPLE_INFO.keys())
MG_SAMPLE_NAMES = ["R-MG", "S-MG", "M-MG"] # Define MG samples for specific analysis
USE_NOVEL_GENE = False
# Novel Gene File Path


# In[ ]:


# Use absolute path or ensure relative path is correct from script execution location
NOVEL_GENE_FILE = "/data02/sunxuebo/project/scrnaseq/no-mammal/mammalnewgene/mammalian_new_genes_final.csv"


# In[ ]:


# Check if novel gene file exists
if not os.path.exists(NOVEL_GENE_FILE):
     warnings.warn(f"Novel gene file not found at: {NOVEL_GENE_FILE}. Skipping novel gene analysis.")
     NOVEL_GENE_FILE = None # Set to None if not found

# Output directory for stage-specific results
OUTPUT_BASE_DIR = "stage_specific_hvg_analysis"
os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)

# --- Helper Functions ---

def get_unique_stages(sample_info_dict):
    """Reads all AnnData objects to find unique stages."""
    unique_stages = set()
    print("Identifying unique stages across all samples...")
    for sample, path in sample_info_dict.items():
        try:
            if os.path.exists(path):
                adata = sc.read(path, cache=True)
                if 'stage' in adata.obs.columns:
                    unique_stages.update(adata.obs['stage'].unique())
                else:
                    print(f"Warning: 'stage' column not found in {sample}")
            else:
                print(f"Warning: File not found for {sample} at {path}")
        except Exception as e:
            print(f"Error reading {sample} from {path}: {e}")
    print(f"Found stages: {list(unique_stages)}")
    return list(unique_stages)

def load_hvgs_for_stage(sample_name, file_path, stage, n_top_genes):
    """Loads data, subsets by stage, finds HVGs."""
    print(f"  Processing {sample_name} for stage '{stage}'...")
    try:
        adata_full = sc.read(file_path, cache=True)

        if 'stage' not in adata_full.obs.columns:
            print(f"  Warning: 'stage' column missing in {sample_name}. Skipping stage subsetting.")
            adata_stage = adata_full # Process full data if no stage info
        else:
             # Check if stage exists in this sample
            if stage not in adata_full.obs['stage'].unique():
                 print(f"  Stage '{stage}' not found in {sample_name}. Skipping sample for this stage.")
                 return None, False # Return None HVGs, Indicate invalid

            adata_stage = adata_full[adata_full.obs['stage'] == stage].copy()

        if adata_stage.n_obs == 0:
            print(f"  No cells found for stage '{stage}' in {sample_name}. Skipping.")
            return None, False # Return None HVGs, Indicate invalid

        print(f"  Found {adata_stage.n_obs} cells for stage '{stage}' in {sample_name}.")

        # Normalize and log-transform (necessary for HVG finding)
        adata_stage.X=adata_stage.layers['counts']
        sc.pp.normalize_total(adata_stage)
        sc.pp.log1p(adata_stage)

        # Find HVGs on the stage-specific subset
        sc.pp.highly_variable_genes(adata_stage, n_top_genes=n_top_genes) # Using seurat_v3 is common
        hvgs = adata_stage.var[adata_stage.var['highly_variable']].index.tolist()
        print(f"  {sample_name} (Stage {stage}): {len(hvgs)} HVGs found.")
        return set(hvgs), True # Return HVG set, Indicate valid

    except FileNotFoundError:
        print(f"  Error: File not found for {sample_name} at {file_path}")
        return None, False
    except Exception as e:
        print(f"  Error processing {sample_name} for stage '{stage}': {e}")
        return None, False

def calculate_intersections_for_plot(hvgs_dict, sample_names, mg_sample_names):
    """Calculates specific intersections needed for the UpSet-style plot."""
    intersections = {}
    valid_sample_names = list(hvgs_dict.keys()) # Samples actually present in this stage

    # 1. All pairwise intersections among valid samples
    for comb in combinations(valid_sample_names, 2):
        intersect_set = hvgs_dict[comb[0]] & hvgs_dict[comb[1]]
        intersections[comb] = len(intersect_set)

    # 2. Three MG sample intersection (if all are present)
    valid_mg_samples = [s for s in mg_sample_names if s in valid_sample_names]
    if len(valid_mg_samples) == 3:
        mg_comb = tuple(sorted(valid_mg_samples)) # Ensure consistent order
        intersect_set = set.intersection(*(hvgs_dict[sample] for sample in mg_comb))
        intersections[mg_comb] = len(intersect_set)
    elif len(valid_mg_samples) > 0 :
         print(f"  Note: Not all 3 MG samples ({mg_sample_names}) present/valid in this stage. Skipping 3-way MG intersection plot point.")

    # 3. Intersection of all valid samples for this stage (if more than 1)
    if len(valid_sample_names) > 1:
         all_comb = tuple(sorted(valid_sample_names))
         intersect_set = set.intersection(*(hvgs_dict[sample] for sample in all_comb))
         intersections[all_comb] = len(intersect_set)

    # Also include single sets for context in the plot if needed (optional, original plot didn't show singles)
    # for sample in valid_sample_names:
    #     intersections[(sample,)] = len(hvgs_dict[sample])

    print(f"  Calculated {len(intersections)} intersection counts for plotting.")
    return intersections

def plot_upset_style(intersections, all_possible_sample_names, output_dir, filename_prefix=""):
    """Plots the UpSet-style visualization for calculated intersections."""
    if not intersections:
        print("  No intersections to plot.")
        return

    # Prepare data - sort by count descending
    combinations_list = list(intersections.keys())
    counts = list(intersections.values())
    # Filter out zero counts before sorting (important!)
    non_zero_indices = [i for i, count in enumerate(counts) if count > 0]
    if not non_zero_indices:
        print("  No non-zero intersections to plot.")
        return

    combinations_list = [combinations_list[i] for i in non_zero_indices]
    counts = [counts[i] for i in non_zero_indices]

    sorted_indices = sorted(range(len(counts)), key=lambda i: counts[i], reverse=True)
    combinations_sorted = [combinations_list[i] for i in sorted_indices]
    counts_sorted = [counts[i] for i in sorted_indices]

    # Determine the samples actually involved in the plotted intersections
    plotted_samples = sorted(list(set(s for comb in combinations_sorted for s in comb)))
    # Use all_possible_sample_names to maintain consistent y-axis order if desired
    # y_axis_samples = all_possible_sample_names
    y_axis_samples = plotted_samples # Or just plot the relevant ones

    # Create figure
    fig = plt.figure(figsize=(max(15, len(combinations_sorted) * 0.8), 8)) # Adjust size dynamically
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[3, 1], hspace=0.05)

    # Bar chart (top)
    ax_bar = fig.add_subplot(grid[0, 0])
    bar_positions = np.arange(len(combinations_sorted))
    bar_width = 0.6
    cmap = plt.cm.get_cmap("viridis", len(combinations_sorted))
    colors = [cmap(i) for i in range(len(combinations_sorted))]
    bars = ax_bar.bar(bar_positions, counts_sorted, width=bar_width, color=colors)
    ax_bar.set_xticks(bar_positions)
    ax_bar.set_xticklabels([]) # Hide x-tick labels
    ax_bar.set_ylabel("Number of Shared HVGs")
    ax_bar.set_title(f"{filename_prefix}Intersection of HVGs ({len(y_axis_samples)} Samples)", fontsize=14)
    ax_bar.set_xlim(-0.5, len(combinations_sorted) - 0.5)
    ax_bar.tick_params(axis='x', length=0) # Hide x-tick marks

    for bar, count in zip(bars, counts_sorted):
        ax_bar.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(5, 0.01 * max(counts_sorted)), # Adjust text offset
                    str(count), ha='center', va='bottom', fontsize=9)

    # Intersection matrix (bottom)
    ax_matrix = fig.add_subplot(grid[1, 0], sharex=ax_bar) # Share x-axis

    # Set alternating background colors for rows
    for i in range(len(y_axis_samples)):
        ax_matrix.axhspan(i - 0.5, i + 0.5, color='lightgrey' if i % 2 == 0 else 'white', alpha=0.6, zorder=0)

    for i, comb in enumerate(combinations_sorted):
        present_indices = [y_axis_samples.index(s) for s in comb if s in y_axis_samples]
        if not present_indices: continue

        # Plot dots
        ax_matrix.scatter([i] * len(present_indices), present_indices, color='black', s=50, zorder=2)
        # Plot connecting lines
        if len(present_indices) > 1:
            ax_matrix.plot([i, i], [min(present_indices), max(present_indices)], color='black', lw=1.5, zorder=1)

    ax_matrix.set_xticks(bar_positions)
    # Set combination labels if needed (can get crowded)
    # comb_labels = [' & '.join(map(str, c)) for c in combinations_sorted]
    # ax_matrix.set_xticklabels(comb_labels, rotation=90, ha='center', fontsize=8)
    ax_matrix.set_xticklabels([]) # Usually better to hide them and rely on the dots

    ax_matrix.set_yticks(np.arange(len(y_axis_samples)))
    ax_matrix.set_yticklabels(y_axis_samples)
    ax_matrix.set_xlim(-0.5, len(combinations_sorted) - 0.5)
    ax_matrix.set_ylim(-0.5, len(y_axis_samples) - 0.5)
    ax_matrix.invert_yaxis() # Place top sample at the top

    plt.tight_layout(rect=[0, 0, 1, 0.97]) # Adjust layout slightly for title
    png_path = os.path.join(output_dir, f"{filename_prefix}HVGs_intersections_upset.png")
    pdf_path = os.path.join(output_dir, f"{filename_prefix}HVGs_intersections_upset.pdf")
    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    plt.savefig(pdf_path, dpi=400, bbox_inches='tight')
    plt.close(fig)
    print(f"  Upset-style plot saved to {png_path} and {pdf_path}")

def calculate_and_export_all_intersections(hvgs_dict, output_dir, filename_prefix=""):
    """Calculates all pairwise and higher-order intersections and saves to CSV."""
    intersections_genes = {}
    valid_sample_names = list(hvgs_dict.keys())

    print("  Calculating all intersections...")
    # Calculate intersections for all combination sizes (2, 3, ..., N)
    for k in range(2, len(valid_sample_names) + 1):
        for comb in combinations(valid_sample_names, k):
            intersect_set = set.intersection(*(hvgs_dict[sample] for sample in comb))
            if intersect_set: # Only store non-empty intersections
                intersection_name = " & ".join(sorted(list(comb)))
                intersections_genes[intersection_name] = sorted(list(intersect_set))
                print(f"    {intersection_name}: {len(intersect_set)} genes")

    if not intersections_genes:
        print("  No intersections found.")
        return None # Indicate no file was saved

    # Export to CSV
    # Pad lists to the same length for DataFrame creation
    max_length = max(len(genes) for genes in intersections_genes.values()) if intersections_genes else 0
    export_data = {name: genes + [""] * (max_length - len(genes)) for name, genes in intersections_genes.items()}
    df = pd.DataFrame(export_data)

    # Sort columns alphabetically for consistency
    df = df[sorted(df.columns)]

    output_csv = os.path.join(output_dir, f"{filename_prefix}HVGs_all_intersections.csv")
    df.to_csv(output_csv, index=False)
    print(f"  All intersection gene lists saved to {output_csv}")
    return output_csv # Return the path to the saved file

def calculate_and_export_unique_hvgs(hvgs_dict, output_dir, filename_prefix=""):
    """Calculates unique HVGs for each sample and saves to CSV."""
    unique_hvgs_dict = {}
    all_unique_data = []
    valid_sample_names = list(hvgs_dict.keys())
    print("  Calculating unique HVGs...")

    for i, name in enumerate(valid_sample_names):
        other_sets = [hvgs_dict[valid_sample_names[j]] for j in range(len(valid_sample_names)) if j != i]
        # Handle case where there are no 'other' sets (only one sample)
        if other_sets:
            union_of_others = set.union(*other_sets)
            unique_hvgs = hvgs_dict[name] - union_of_others
        else:
            unique_hvgs = hvgs_dict[name] # If only one sample, all its HVGs are unique

        unique_hvgs_dict[name] = unique_hvgs
        print(f"    Unique HVGs for {name}: {len(unique_hvgs)}")
        for hvg in sorted(list(unique_hvgs)):
            all_unique_data.append([name, hvg])

    if not all_unique_data:
        print("  No unique HVGs found.")
        return None, None # Indicate no file was saved and return empty dict

    df = pd.DataFrame(all_unique_data, columns=["Sample", "HVG"])
    output_csv = os.path.join(output_dir, f"{filename_prefix}HVGs_unique.csv")
    df.to_csv(output_csv, index=False)
    print(f"  Unique HVGs saved to {output_csv}")
    return output_csv, unique_hvgs_dict # Return path and the dictionary

def perform_mg_specific_analysis(hvgs_dict, mg_sample_names, output_dir):
    """Performs Venn, Bar chart, and intersection saving specific to MG samples."""
    valid_mg_samples = [s for s in mg_sample_names if s in hvgs_dict]

    if len(valid_mg_samples) < 2:
        print(f"  Skipping MG-specific analysis: Found {len(valid_mg_samples)}/{len(mg_sample_names)} MG samples valid for this stage.")
        return

    print(f"  Performing MG-specific analysis for: {valid_mg_samples}")
    mg_hvgs_dict = {s: hvgs_dict[s] for s in valid_mg_samples}

    # 1. Venn Diagram (only works well for 2 or 3 sets)
    if len(valid_mg_samples) in [2, 3]:
        try:
            plt.figure(figsize=(8, 8))
            venn(mg_hvgs_dict)
            plt.title(f"Venn Diagram of MG HVGs ({' & '.join(valid_mg_samples)})")
            venn_png_path = os.path.join(output_dir, "MG_HVGs_Venn.png")
            venn_pdf_path = os.path.join(output_dir, "MG_HVGs_Venn.pdf")
            plt.savefig(venn_png_path, dpi=300, bbox_inches='tight')
            plt.savefig(venn_pdf_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"    MG Venn diagram saved to {venn_png_path} and {venn_pdf_path}")
        except Exception as e:
            print(f"    Error generating Venn diagram: {e}") # Venn lib can sometimes fail

    # 2. Bar Chart of Pairwise and (if 3 samples) Triple Intersections
    shared_counts = {}
    shared_genes = {} # Store actual genes for saving later

    # Pairwise
    for s1, s2 in combinations(valid_mg_samples, 2):
        shared = mg_hvgs_dict[s1] & mg_hvgs_dict[s2]
        label = f"{s1} & {s2}"
        shared_counts[label] = len(shared)
        shared_genes[label] = shared
        print(f"      Shared {label}: {len(shared)}")

    # Triple (if 3 samples)
    if len(valid_mg_samples) == 3:
        shared_all = set.intersection(*mg_hvgs_dict.values())
        label_all = f"All ({' & '.join(valid_mg_samples)})"
        shared_counts[label_all] = len(shared_all)
        shared_genes[label_all] = shared_all
        print(f"      Shared All 3 MG: {len(shared_all)}")

    if shared_counts:
        plt.figure(figsize=(8, 6))
        categories = list(shared_counts.keys())
        values = list(shared_counts.values())
        # Use distinct colors
        colors = plt.cm.get_cmap('Pastel1', len(categories)).colors
        plt.bar(categories, values, color=colors, width=0.5)
        plt.xlabel("MG Sample Groups", fontsize=12)
        plt.ylabel("Number of Shared HVGs", fontsize=12)
        plt.title("Shared HVGs Among MG Samples", fontsize=14)
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        bar_png_path = os.path.join(output_dir, "MG_HVGs_shared_barplot.png")
        bar_pdf_path = os.path.join(output_dir, "MG_HVGs_shared_barplot.pdf")
        plt.savefig(bar_png_path, dpi=300, bbox_inches='tight')
        plt.savefig(bar_pdf_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    MG shared HVGs bar plot saved to {bar_png_path} and {bar_pdf_path}")

        # 3. Save MG intersection genes to CSV
        mg_intersect_data = []
        for group, genes in shared_genes.items():
            for gene in sorted(list(genes)):
                 mg_intersect_data.append([group, gene])

        if mg_intersect_data:
             df_mg_intersect = pd.DataFrame(mg_intersect_data, columns=["Intersection Group", "Gene"])
             mg_csv_path = os.path.join(output_dir, "MG_HVGs_intersections.csv")
             df_mg_intersect.to_csv(mg_csv_path, index=False)
             print(f"    MG intersection genes saved to {mg_csv_path}")


def intersect_with_novel_genes(hvg_sets_dict, novel_genes_set, set_type_name, output_dir, filename_prefix=""):
    """
    Intersects sets of HVGs (e.g., intersections, unique) with a set of novel genes.

    Args:
        hvg_sets_dict (dict): Dict where keys are group names (e.g., "SampleA & SampleB", "SampleC Unique")
                              and values are sets or lists of gene names.
        novel_genes_set (set): Set of novel gene names.
        set_type_name (str): Type of HVG set (e.g., "intersection", "unique") for filename.
        output_dir (str): Directory to save the results.
        filename_prefix (str): Optional prefix for the output filename.
    """
    if not novel_genes_set:
        print(f"  Skipping novel gene intersection for {set_type_name} sets: Novel gene list not available.")
        return

    print(f"  Intersecting {set_type_name} HVGs with novel genes...")
    novel_intersections = {}
    for group_name, gene_set in hvg_sets_dict.items():
        # Ensure we are working with a set
        if isinstance(gene_set, list):
            gene_set = set(gene_set)
        elif not isinstance(gene_set, set):
             print(f"    Warning: Unexpected data type for {group_name}. Skipping.")
             continue

        novel_found = list(gene_set.intersection(novel_genes_set))
        if novel_found:
            novel_intersections[group_name] = sorted(novel_found)
            print(f"    {group_name}: Found {len(novel_found)} novel genes.")
        # else:
            # print(f"    {group_name}: No novel genes found.")


    if not novel_intersections:
        print(f"  No novel genes found in any {set_type_name} set.")
        return

    # Export to CSV
    max_length = max(len(genes) for genes in novel_intersections.values()) if novel_intersections else 0
    export_data = {name: genes + [""] * (max_length - len(genes)) for name, genes in novel_intersections.items()}
    df = pd.DataFrame(export_data)

     # Sort columns alphabetically for consistency
    df = df[sorted(df.columns)]

    output_csv = os.path.join(output_dir, f"{filename_prefix}HVGs_{set_type_name}_novel_gene_overlap.csv")
    df.to_csv(output_csv, index=False)
    print(f"  Novel gene overlap for {set_type_name} sets saved to {output_csv}")


# --- Main Execution Logic ---

def main():
    # 1. Load Novel Genes (once)
    novel_genes_set = set()
    if NOVEL_GENE_FILE and USE_NOVEL_GENE:
        try:
            novel_df = pd.read_csv(NOVEL_GENE_FILE)
            # Assuming the column with relevant genes is 'mouse' - **ADJUST IF NEEDED**
            if 'mouse' in novel_df.columns:
                 novel_genes_set = set(novel_df["mouse"].dropna().astype(str).tolist())
                 print(f"Loaded {len(novel_genes_set)} unique novel genes for comparison.")
            else:
                 print(f"Warning: Column 'mouse' not found in {NOVEL_GENE_FILE}. Cannot load novel genes.")
        except Exception as e:
            print(f"Error loading novel gene file {NOVEL_GENE_FILE}: {e}")
    else:
         print("Novel gene file not specified or not found. Skipping novel gene analysis.")


    # 2. Identify Unique Stages
    stages = get_unique_stages(SAMPLE_INFO)
    if not stages:
        print("No stages identified or error reading data. Exiting.")
        return

    # 3. Process Each Stage
    for stage in stages:
        print(f"\n--- Processing Stage: {stage} ---")
        stage_output_dir = os.path.join(OUTPUT_BASE_DIR, str(stage))
        os.makedirs(stage_output_dir, exist_ok=True)

        # 3.1 Load HVGs for this stage
        stage_hvgs_dict = {}
        stage_valid_samples = []
        for sample_name in ALL_SAMPLE_NAMES:
            file_path = SAMPLE_INFO[sample_name]
            hvgs, is_valid = load_hvgs_for_stage(sample_name, file_path, stage, N_TOP_GENES)
            if is_valid:
                stage_hvgs_dict[sample_name] = hvgs
                stage_valid_samples.append(sample_name)

        if len(stage_valid_samples) < 2:
             print(f"Stage '{stage}': Found {len(stage_valid_samples)} valid samples. Need at least 2 for comparative analysis. Skipping stage.")
             continue

        print(f"\nStage '{stage}': Analysis based on valid samples: {stage_valid_samples}")

        # 3.2 Calculate Intersections for UpSet-style plot
        intersections_for_plot = calculate_intersections_for_plot(stage_hvgs_dict, stage_valid_samples, MG_SAMPLE_NAMES)

        # 3.3 Plot UpSet-style plot
        plot_upset_style(intersections_for_plot, stage_valid_samples, stage_output_dir, filename_prefix=f"Stage_{stage}_")

        # 3.4 Calculate and Export *All* Intersections (pairwise, triplets, etc.)
        all_intersections_csv = calculate_and_export_all_intersections(stage_hvgs_dict, stage_output_dir, filename_prefix=f"Stage_{stage}_")

        # 3.5 Calculate and Export Unique HVGs
        unique_hvgs_csv, unique_hvgs_dict = calculate_and_export_unique_hvgs(stage_hvgs_dict, stage_output_dir, filename_prefix=f"Stage_{stage}_")

        # 3.6 Perform MG-Specific Analysis (if applicable)
        perform_mg_specific_analysis(stage_hvgs_dict, MG_SAMPLE_NAMES, stage_output_dir) # Output filenames inside function are relative to stage_output_dir

        # 3.7 Intersect with Novel Genes (if available)
        if novel_genes_set:
            # Intersect ALL intersection sets with novel genes
            if all_intersections_csv: # Check if the file was actually created
                 try:
                     # Read the saved intersection data back (or pass the dictionary directly if refactored)
                      intersections_df_for_novel = pd.read_csv(all_intersections_csv)
                      intersection_sets_for_novel = {col: set(intersections_df_for_novel[col].dropna().astype(str))
                                                    for col in intersections_df_for_novel.columns}
                      intersect_with_novel_genes(intersection_sets_for_novel, novel_genes_set, "intersection", stage_output_dir, filename_prefix=f"Stage_{stage}_")
                 except Exception as e:
                     print(f"  Error processing intersections CSV for novel gene overlap: {e}")


            # Intersect UNIQUE gene sets with novel genes
            if unique_hvgs_dict: # Check if unique HVGs were found
                 intersect_with_novel_genes(unique_hvgs_dict, novel_genes_set, "unique", stage_output_dir, filename_prefix=f"Stage_{stage}_")


    print("\n--- Analysis Complete ---")

if __name__ == "__main__":
    # Suppress excessive warnings from dependencies if desired
    # warnings.filterwarnings('ignore')
    main()


# In[ ]:


# 1. 包含各个 stage 子目录的基础路径
#    假设你的脚本和 'stage_specific_hvg_analysis' 目录在同一级
#    或者提供绝对路径 '/path/to/your/stage_specific_hvg_analysis'
BASE_OUTPUT_DIR = "./stage_specific_hvg_analysis"

# 2. 需要处理的 stage 名称列表
STAGES_TO_PROCESS = ['stage0', 'stage1', 'stage2', 'stage3', 'stage4', 'ALL']
#    注意：'ALL' 可能是你对所有阶段合并分析的结果目录名，请根据实际情况修改

# 4. 输出整合后的CSV文件名
CONSOLIDATED_OUTPUT_FILE = "consolidated_target_novel_genes.csv"

# --- 脚本执行 ---

# 用于存储最终结果的列表，每个元素是一个包含 [stage, intersection, gene] 的列表
consolidated_data = []

print("开始处理各个 stage 的 novel gene intersection 文件...")

# 遍历每个指定的 stage
for stage_name in STAGES_TO_PROCESS:
    print(f"\n--- 正在处理 Stage: {stage_name} ---")

    # 构建当前 stage 的子目录路径
    stage_dir = os.path.join(BASE_OUTPUT_DIR, stage_name)

    # 构建目标 CSV 文件的完整路径
    # 文件名格式根据你的描述是 'Stage_{stage_name}_HVGs_intersection_novel_gene_overlap.csv'
    if USE_NOVEL_GENE:
        csv_filename = f'Stage_{stage_name}_HVGs_intersection_novel_gene_overlap.csv'
    else:
        csv_filename = f'Stage_{stage_name}_HVGs_all_intersections.csv'
    csv_filepath = os.path.join(stage_dir, csv_filename)

    # 检查目标 CSV 文件是否存在
    if not os.path.exists(csv_filepath):
        warnings.warn(f"文件未找到，跳过: {csv_filepath}")
        continue # 跳到下一个 stage

    # 读取 CSV 文件
    try:
        df = pd.read_csv(csv_filepath)
        print(f"  成功读取文件: {csv_filename}")

        # 遍历需要提取的目标交集列名
        for intersection_col in TARGET_INTERSECTIONS:
            # 检查目标列是否存在于当前 DataFrame 中
            if intersection_col in df.columns:
                print(f"    找到目标列: '{intersection_col}'")
                # 提取该列的基因名
                # 使用 .dropna() 去除因为列长不一致而填充的 NaN/空值
                # 使用 .astype(str) 确保基因名是字符串
                # 使用 .unique() 可选，如果原始文件可能有重复基因且你只想保留唯一值
                genes_in_intersection = df[intersection_col].dropna().astype(str).tolist()

                if genes_in_intersection:
                    print(f"      提取到 {len(genes_in_intersection)} 个基因.")
                    # 将每个基因添加到结果列表中
                    for gene in genes_in_intersection:
                        # 确保基因名不是空字符串（虽然dropna通常能处理，加一层保险）
                        if gene.strip():
                            consolidated_data.append([stage_name, intersection_col, gene.strip()])
                else:
                    print(f"      列 '{intersection_col}' 中没有找到有效基因名。")

            else:
                # 如果目标列不存在，发出警告
                warnings.warn(f"在文件 {csv_filename} 中未找到列: '{intersection_col}'")

    except Exception as e:
        # 处理读取文件时可能发生的其他错误
        warnings.warn(f"读取文件时发生错误 {csv_filepath}: {e}")
        continue # 跳到下一个 stage

# 将收集到的数据转换为 DataFrame
if consolidated_data:
    consolidated_df = pd.DataFrame(consolidated_data, columns=['stage', 'intersection', 'gene'])
    print(f"\n数据整合完成，共收集到 {len(consolidated_df)} 条记录。")

    # 保存 DataFrame 到 CSV 文件
    try:
        consolidated_df.to_csv(CONSOLIDATED_OUTPUT_FILE, index=False)
        print(f"结果已成功保存到: {CONSOLIDATED_OUTPUT_FILE}")
    except Exception as e:
        print(f"保存最终CSV文件时出错: {e}")

else:
    print("\n未能收集到任何数据，未生成输出文件。请检查输入目录、文件名和列名是否正确。")

print("\n脚本执行完毕。")

