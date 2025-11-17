"""
Modular Scanpy workflow for the bone marrow / PBMC dataset.

Usage:
    python bonemarrow_stage2.py --input-path bone_marrow.h5ad --output-h5ad processed.h5ad

Key steps:
1. Load AnnData with error handling.
2. Run quality control, filtering, normalization, and log-transform.
3. Select highly variable genes, scale, perform PCA/UMAP, and cluster.
4. Detect marker genes and annotate clusters using reference `Cell.group`.
5. Optionally visualize intermediate results or write processed data to disk.
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, Optional

import anndata as ad
import pandas as pd
import scanpy as sc


def load_anndata(path: Path) -> ad.AnnData:
    """Load an AnnData object, raising a readable error if unavailable."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    logging.info("Reading AnnData from %s", path)
    try:
        adata = sc.read_h5ad(path)
    except OSError as exc:
        raise OSError(f"Failed to read {path}: {exc}") from exc
    logging.info("Loaded dataset with %d cells and %d genes", adata.n_obs, adata.n_vars)
    return adata


def add_qc_metrics(adata: ad.AnnData) -> None:
    """Attach QC metrics (mitochondrial content, counts) to the object."""
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)


def filter_cells(
    adata: ad.AnnData,
    min_genes: int = 300,
    max_genes: int = 4000,
    max_pct_mt: float = 10.0,
) -> ad.AnnData:
    """Filter out likely empty droplets, doublets, or stressed cells."""
    mask = (
        (adata.obs["n_genes_by_counts"] > min_genes)
        & (adata.obs["n_genes_by_counts"] < max_genes)
        & (adata.obs["pct_counts_mt"] < max_pct_mt)
    )
    removed = (~mask).sum()
    logging.info("Filtering removed %d cells (%.2f%%)", removed, removed / adata.n_obs * 100)
    return adata[mask].copy()


def normalize_and_log(adata: ad.AnnData) -> None:
    """Normalize library size and log-transform counts."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # preserve log-normalized values for marker calling


def select_hvgs(adata: ad.AnnData) -> ad.AnnData:
    """Restrict data frame to highly variable genes."""
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, flavor="seurat_v3"
    )
    hvgs = adata[:, adata.var["highly_variable"]].copy()
    logging.info("Retained %d highly variable genes", hvgs.n_vars)
    return hvgs


def scale_and_embed(adata: ad.AnnData, n_pcs: int = 40, n_neighbors: int = 10) -> None:
    """Scale genes, compute PCA, neighbors, and UMAP."""
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)


def cluster_and_rank(adata: ad.AnnData, resolution: float = 0.5) -> None:
    """Perform Leiden clustering and compute marker genes."""
    sc.tl.leiden(adata, resolution=resolution)
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")


def annotate_cell_types(
    adata: ad.AnnData, reference_key: str = "Cell.group", cluster_key: str = "leiden"
) -> Dict[str, str]:
    """Map clusters to cell types via majority vote using reference annotations."""
    if reference_key not in adata.obs.columns:
        raise KeyError(f"{reference_key} not found in adata.obs")
    mapping = (
        adata.obs.groupby(cluster_key)[reference_key]
        .agg(lambda x: x.value_counts().idxmax())
        .to_dict()
    )
    adata.obs["cell_type"] = adata.obs[cluster_key].map(mapping).astype("category")
    return mapping


def summarize_cell_types(adata: ad.AnnData, key: str = "cell_type") -> pd.DataFrame:
    """Return counts and percentages for a given annotation key."""
    counts = adata.obs[key].value_counts().sort_index()
    df = pd.DataFrame({"n_cells": counts, "percent": counts / counts.sum() * 100})
    logging.info("Cell-type summary:\n%s", df)
    return df


def plot_results(
    adata: ad.AnnData,
    plots_dir: Optional[Path] = None,
    show: bool = False,
) -> None:
    """Optionally render or save diagnostic plots."""
    if plots_dir:
        plots_dir.mkdir(parents=True, exist_ok=True)

    def _save(fig, name: str) -> None:
        if fig is not None and plots_dir:
            fig.savefig(plots_dir / f"{name}.png", dpi=150, bbox_inches="tight")

    fig = sc.pl.pca_variance_ratio(adata, log=True, show=show, return_fig=True)
    _save(fig, "pca_variance_ratio")

    fig = sc.pl.umap(
        adata, color=["n_genes_by_counts", "pct_counts_mt"], show=show, return_fig=True
    )
    _save(fig, "umap_qc")

    fig = sc.pl.umap(adata, color=["leiden"], show=show, return_fig=True)
    _save(fig, "umap_leiden")

    if "cell_type" in adata.obs:
        fig = sc.pl.umap(adata, color=["cell_type"], legend_loc="on data", show=show, return_fig=True)
        _save(fig, "umap_cell_types")

    fig = sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False, show=show, return_fig=True)
    _save(fig, "marker_genes")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the Scanpy pipeline on a bone marrow/PBMC dataset.")
    parser.add_argument("--input-path", type=Path, default=Path("bone_marrow.h5ad"), help="Path to the input .h5ad file.")
    parser.add_argument("--output-h5ad", type=Path, help="Optional path to write the processed AnnData object.")
    parser.add_argument("--plots-dir", type=Path, help="Directory to save figures (PNG).")
    parser.add_argument("--show-plots", action="store_true", help="Display plots interactively.")
    parser.add_argument("--reference-label", default="Cell.group", help="obs column used for majority-vote annotation.")
    parser.add_argument("--min-genes", type=int, default=300, help="Minimum genes per cell for filtering.")
    parser.add_argument("--max-genes", type=int, default=4000, help="Maximum genes per cell for filtering.")
    parser.add_argument("--max-pct-mt", type=float, default=10.0, help="Maximum mitochondrial percentage per cell.")
    parser.add_argument("--resolution", type=float, default=0.5, help="Leiden resolution parameter.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    logging.basicConfig(level=logging.INFO, format="[%(levelname)s] %(message)s")
    try:
        adata = load_anndata(args.input_path)
        add_qc_metrics(adata)
        adata = filter_cells(adata, args.min_genes, args.max_genes, args.max_pct_mt)
        normalize_and_log(adata)
        hvgs = select_hvgs(adata)
        scale_and_embed(hvgs)
        cluster_and_rank(hvgs, resolution=args.resolution)
        mapping = annotate_cell_types(hvgs, reference_key=args.reference_label)
        logging.info("Cluster-to-cell-type mapping: %s", mapping)
        summary = summarize_cell_types(hvgs)
        print(summary)  # explicit stdout for notebooks/CLI capture
        plot_results(hvgs, plots_dir=args.plots_dir, show=args.show_plots)
        if args.output_h5ad:
            hvgs.write(args.output_h5ad)
            logging.info("Wrote processed AnnData to %s", args.output_h5ad)
    except Exception as exc:  # pylint: disable=broad-except
        logging.error("Pipeline failed: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main()
