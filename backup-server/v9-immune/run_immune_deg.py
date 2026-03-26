import os
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# ----------------------------
# Config
# ----------------------------
datasetlist = ['M-MG', 'R-MG', 'S-MG', 'R-AG', 'R-CG', 'S-AG']
group_names = ['newcelltype', 'celltype']
BASE_OUT = './deg-result-imm'

# If you want tighter control over scanpy figure defaults:
# sc.set_figure_params(dpi=120, dpi_save=300)


def ensure_dirs(group_name: str):
    plot_dir = os.path.join(BASE_OUT, group_name, 'plots')
    table_dir = os.path.join(BASE_OUT, group_name, 'table')
    os.makedirs(plot_dir, exist_ok=True)
    os.makedirs(table_dir, exist_ok=True)
    return plot_dir, table_dir


# ----------------------------
# UMAP (PDF, rasterized)
# ----------------------------
def save_umap_pdf_rasterized(adata, color: str, out_pdf: str, dpi: int = 300):
    """Save a UMAP as PDF with rasterized points (keeps text/vector, rasterizes scatter)."""
    # Assumes UMAP already exists (no fallback/compute here)
    if 'X_umap' not in adata.obsm_keys():
        raise KeyError("UMAP not found: adata.obsm['X_umap'] is missing")

    fig = sc.pl.umap(adata, color=color, show=False, return_fig=True)
    # Rasterize any point collections
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    fig.savefig(out_pdf, format='pdf', dpi=dpi, bbox_inches='tight')
    plt.close(fig)


# ----------------------------
# Barplot (stage x celltypename)
# ----------------------------
def do_cell_barplot(ea, dataset: str, celltypename: str, out_png: str, out_pdf: str):
    """Stacked horizontal barplot of percentages over stage for a given celltypename."""
    if 'stage' not in ea.obs.columns:
        raise KeyError("Missing required obs column: 'stage'")
    if celltypename not in ea.obs.columns:
        raise KeyError(f"Missing required obs column: '{celltypename}'")

    # Ensure categorical
    if not pd.api.types.is_categorical_dtype(ea.obs[celltypename]):
        ea.obs[celltypename] = ea.obs[celltypename].astype('category')

    # 1) counts & percentages
    celltype_counts = ea.obs.groupby(['stage', celltypename]).size().unstack(fill_value=0)
    celltype_percentages = celltype_counts.div(celltype_counts.sum(axis=1), axis=0) * 100

    # 2) color map from AnnData uns
    categories = list(ea.obs[celltypename].cat.categories)
    color_map = None
    uns_key = f'{celltypename}_colors'
    if uns_key in ea.uns:
        color_list = ea.uns[uns_key]
        if len(color_list) == len(categories):
            color_map = dict(zip(categories, color_list))

    # If no valid color map, fallback to matplotlib default cycle (still deterministic)
    if color_map is None:
        default_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', [])
        if not default_cycle:
            default_cycle = ['#808080']
        color_map = {ct: default_cycle[i % len(default_cycle)] for i, ct in enumerate(categories)}

    # 3) prepare plot order
    celltype_percentages = celltype_percentages.sort_index(ascending=False)
    times = celltype_percentages.index

    cell_types_to_plot = [ct for ct in categories if ct in celltype_percentages.columns]
    bottom = pd.Series([0] * len(times), index=times)

    # 4) plot
    fig, ax = plt.subplots(figsize=(8, 4))
    bar_width = 0.8

    for cell_type in cell_types_to_plot:
        percentages = celltype_percentages[cell_type]
        ax.barh(
            times,
            percentages,
            left=bottom,
            color=color_map.get(cell_type, '#808080'),
            label=cell_type,
            height=bar_width,
        )
        bottom += percentages

    # 5) format & save
    ax.set_title(f'Percentage of Each Celltype Over Stage in {dataset}', fontsize=14)
    ax.set_xlabel('Percentage', fontsize=12)
    ax.set_ylabel('Stage', fontsize=12)
    ax.legend(title=celltypename, bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches='tight')
    fig.savefig(out_pdf, bbox_inches='tight')
    plt.close(fig)


# ----------------------------
# DEG per group_name
# ----------------------------
def do_deg_one_group(adata, dataset: str, group_name: str):
    plot_dir, table_dir = ensure_dirs(group_name)

    if group_name not in adata.obs.columns:
        print(f"[SKIP] {dataset}: obs does not contain '{group_name}'")
        return

    # DEG (Wilcoxon) on current adata.X (already normalized)
    sc.tl.rank_genes_groups(adata, groupby=group_name, method='wilcoxon')

    # Dotplot -> custom path
    dot_png = os.path.join(plot_dir, f'{dataset}_{group_name}_rank_genes_dotplot.png')
    dp = sc.pl.rank_genes_groups_dotplot(
        adata,
        groupby=group_name,
        n_genes=5,
        min_logfoldchange=0.25,
        show=False,
        return_fig=True,
    )
    dp.savefig(dot_png, dpi=300, bbox_inches='tight')
    plt.close('all')

    # DEG table
    deg_csv = os.path.join(table_dir, f'{dataset}_{group_name}_ranked_genes_for_annotation.csv')
    deg_df = sc.get.rank_genes_groups_df(
        adata,
        group=None,
        pval_cutoff=0.01,
        log2fc_min=0.1,
    )
    deg_df.to_csv(deg_csv, index=False)

    # UMAP (PDF, rasterized)
    umap_pdf = os.path.join(plot_dir, f'{dataset}_umap_color_{group_name}.pdf')
    save_umap_pdf_rasterized(adata, color=group_name, out_pdf=umap_pdf)

    # Barplot (png+pdf)
    bar_png = os.path.join(plot_dir, f'{dataset}_{group_name}_bar_plot.png')
    bar_pdf = os.path.join(plot_dir, f'{dataset}_{group_name}_bar_plot.pdf')
    do_cell_barplot(adata, dataset=dataset, celltypename=group_name, out_png=bar_png, out_pdf=bar_pdf)

    print(f"[OK] {dataset} / {group_name} -> {plot_dir} | {table_dir}")


def run_all():
    for dataset in datasetlist:
        adata_path = f'./{dataset}/{dataset}_lum-imm_raw.h5ad'
        if not os.path.exists(adata_path):
            print(f"[SKIP] Missing file: {adata_path}")
            continue

        print(f"[LOAD] {dataset}: {adata_path}")
        adata = sc.read_h5ad(adata_path)

        # No fallback/compute: require precomputed UMAP
        if 'X_umap' not in adata.obsm_keys():
            raise KeyError(f"{dataset}: Missing UMAP (adata.obsm['X_umap'])")

        for group_name in group_names:
            do_deg_one_group(adata, dataset, group_name)


if __name__ == '__main__':
    run_all()

