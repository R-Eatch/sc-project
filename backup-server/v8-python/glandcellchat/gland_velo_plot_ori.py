#!/usr/bin/env python
# coding: utf-8
"""Plot scVelo stream plots with legend on the right margin.

This script assumes the velocity computation has already finished and the
processed AnnData objects were saved by your existing pipeline.

Default inputs (per dataset):
  ./data/{DATASET}_velofinished.h5ad

Outputs:
  ./figure_velo/{DATASET}_velo_rightLABEL.pdf

Notes:
- If you only have *_pp.h5ad (no velocity results), this script will not
  recompute velocity; it will fail because velocity_graph is missing.
"""

import warnings
warnings.filterwarnings("ignore")

import os
import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import scvelo as scv

# --- scVelo Settings ---
scv.logging.print_version()
scv.settings.verbosity = 3
scv.settings.presenter_view = True

# Where to save figures
scv.settings.figdir = "./figure_velo_ori"
os.makedirs(scv.settings.figdir, exist_ok=True)

# --- Global Variables (keep consistent with your pipeline) ---
glandlist = ["MG", "SG", "EG"]

# Which categorical label to color by
# If you used newcelltype/leiden instead, change this.
color_key = "celltype"

# Where to read finished objects from
in_tmpl = "./data/{ds}_velofinished_ori.h5ad"


def rasterize_all_artists(fig: matplotlib.figure.Figure) -> None:
    """Rasterize artists to keep PDF size reasonable for dense plots."""
    for artist in fig.findobj(match=matplotlib.artist.Artist):
        if hasattr(artist, "set_rasterized"):
            try:
                artist.set_rasterized(True)
            except Exception:
                pass


def plot_right_legend(ds: str) -> None:
    print(f"--- Plotting {ds} ---")

    file_path = in_tmpl.format(ds=ds)
    if not os.path.exists(file_path):
        print(f"File not found: {file_path} (skipping)")
        return

    adata = sc.read_h5ad(file_path)

    # Sanity checks
    if "X_umap" not in adata.obsm:
        print(f"X_umap not found in {ds}. Cannot plot velocity stream (skipping)")
        return

    # velocity_graph is required for velocity_embedding_stream
    if "velocity_graph" not in adata.uns:
        # Some scVelo versions store it in .uns or .obsp; check both.
        # If missing, we cannot plot without recomputing.
        has_obsp = any(k.startswith("velocity_graph") for k in getattr(adata, "obsp", {}).keys())
        if not has_obsp:
            print(
                f"velocity_graph not found in {ds}. "
                "This file may not be a *_velofinished.h5ad from your pipeline. Skipping."
            )
            return

    if color_key not in adata.obs:
        # Fallbacks commonly used in your scripts
        for alt in ["newcelltype", "leiden"]:
            if alt in adata.obs:
                print(f"Warning: '{color_key}' not in adata.obs for {ds}; using '{alt}' instead")
                plot_color = alt
                break
        else:
            print(f"No suitable categorical label found in adata.obs for {ds} (skipping)")
            return
    else:
        plot_color = color_key

    ax = scv.pl.velocity_embedding_stream(
        adata,
        basis="umap",
        color=plot_color,
        legend_loc="right margin",  # key: move legend to right sidebar
        show=False,
        title=f"{ds} Velocity",
    )

    if ax is None:
        print(f"Plot failed (ax is None) for {ds} (skipping)")
        return

    fig = ax.get_figure()
    rasterize_all_artists(fig)

    out_pdf = os.path.join(scv.settings.figdir, f"{ds}_velo_rightLABEL_ori.pdf")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight", format="pdf")
    plt.close(fig)

    print(f"Saved: {out_pdf}")


def main():
    for ds in glandlist:
        plot_right_legend(ds)


if __name__ == "__main__":
    main()

