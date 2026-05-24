from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmwrite


# -----------------------------
# 1. Paths
# -----------------------------
base_dir = Path(__file__).resolve().parents[1]

sim_dir = base_dir / "data" / "simulated"
out_base = base_dir / "data" / "decontx_input"
figure_dir = base_dir / "results" / "figures"

out_base.mkdir(parents=True, exist_ok=True)
figure_dir.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = figure_dir
sc.settings.verbosity = 2


# -----------------------------
# 2. Helper
# -----------------------------
def get_counts(adata):
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def make_precluster_labels(adata, condition):
    """
    Generate DecontX input clusters using only contaminated observed counts.
    This avoids using clean ground-truth cluster labels.
    """
    adata_tmp = adata.copy()

    counts = get_counts(adata_tmp)
    adata_tmp.X = counts.copy()

    if "counts" not in adata_tmp.layers:
        adata_tmp.layers["counts"] = adata_tmp.X.copy()

    # Standard Scanpy clustering using contaminated counts
    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)

    sc.pp.highly_variable_genes(
        adata_tmp,
        n_top_genes=2000,
        flavor="seurat"
    )

    sc.pp.pca(
        adata_tmp,
        use_highly_variable=True,
        svd_solver="arpack"
    )

    sc.pp.neighbors(
        adata_tmp,
        n_neighbors=10,
        n_pcs=40
    )

    sc.tl.umap(adata_tmp)

    sc.tl.leiden(
        adata_tmp,
        resolution=0.5,
        key_added="decontx_input_cluster"
    )

    # Save UMAP for checking the pre-clustering used by DecontX
    sc.pl.umap(
        adata_tmp,
        color=["decontx_input_cluster"],
        title=f"DecontX input clusters: {condition}",
        save=f"_decontx_input_clusters_{condition}.png",
        show=False,
    )

    return adata_tmp.obs["decontx_input_cluster"].astype(str).values


# -----------------------------
# 3. Find simulated datasets
# -----------------------------
sim_files = sorted(sim_dir.glob("pbmc10k_contaminated_*.h5ad"))

if len(sim_files) == 0:
    raise FileNotFoundError(f"No simulated h5ad files found in: {sim_dir}")

print("Found simulated files:")
for f in sim_files:
    print(" -", f.name)


# -----------------------------
# 4. Export one folder per condition
# -----------------------------
summary_records = []

for sim_file in sim_files:
    condition = sim_file.stem.replace("pbmc10k_contaminated_", "")
    print(f"\nExporting DecontX input for: {condition}")

    adata = sc.read_h5ad(sim_file)
    counts = get_counts(adata)

    if not sparse.issparse(counts):
        counts = sparse.csr_matrix(counts)

    # Important:
    # AnnData stores cells x genes.
    # DecontX / SingleCellExperiment expects genes x cells.
    counts_genes_by_cells = counts.T.tocoo()

    out_dir = out_base / condition
    out_dir.mkdir(parents=True, exist_ok=True)

    # Make cluster labels from the contaminated data itself
    decontx_clusters = make_precluster_labels(adata, condition)

    # Save Matrix Market format
    mmwrite(out_dir / "counts_genes_by_cells.mtx", counts_genes_by_cells)

    # Save gene names and cell barcodes separately
    pd.DataFrame({"gene": adata.var_names}).to_csv(
        out_dir / "genes.csv",
        index=False
    )

    pd.DataFrame({"barcode": adata.obs_names}).to_csv(
        out_dir / "barcodes.csv",
        index=False
    )

    pd.DataFrame({
        "barcode": adata.obs_names,
        "decontx_input_cluster": decontx_clusters
    }).to_csv(
        out_dir / "clusters.csv",
        index=False
    )

    cluster_counts = pd.Series(decontx_clusters).value_counts().sort_index()

    for cluster, n_cells in cluster_counts.items():
        summary_records.append({
            "condition": condition,
            "cluster": cluster,
            "n_cells": n_cells
        })

    print("Saved:")
    print(" ", out_dir / "counts_genes_by_cells.mtx")
    print(" ", out_dir / "genes.csv")
    print(" ", out_dir / "barcodes.csv")
    print(" ", out_dir / "clusters.csv")


summary_df = pd.DataFrame(summary_records)
summary_path = out_base / "decontx_input_cluster_summary.csv"
summary_df.to_csv(summary_path, index=False)

print("\nSaved cluster summary:")
print(summary_path)

print("\nDone exporting optimized DecontX inputs.")
