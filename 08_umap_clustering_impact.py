from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Paths
# -----------------------------
base_dir = Path(__file__).resolve().parents[1]

clean_path = base_dir / "data" / "processed" / "pbmc10k_clean_baseline.h5ad"
sim_dir = base_dir / "data" / "simulated"
decontx_dir = base_dir / "data" / "decontx_output"
decontx_input_dir = base_dir / "data" / "decontx_input"

figure_dir = base_dir / "results" / "figures"
table_dir = base_dir / "results" / "tables"

figure_dir.mkdir(parents=True, exist_ok=True)
table_dir.mkdir(parents=True, exist_ok=True)

sc.settings.verbosity = 2


# -----------------------------
# Helper functions
# -----------------------------
def as_csr(x):
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def get_counts(adata):
    if "counts" in adata.layers:
        return as_csr(adata.layers["counts"])
    return as_csr(adata.X)


def condition_order_key(condition):
    condition = str(condition)

    if condition.startswith("1gene"):
        ambient_rank = 0
    elif condition.startswith("100genes"):
        ambient_rank = 1
    else:
        ambient_rank = 2

    if "5pct" in condition and "15pct" not in condition:
        contam_rank = 0
    elif "15pct" in condition:
        contam_rank = 1
    elif "30pct" in condition:
        contam_rank = 2
    else:
        contam_rank = 3

    return ambient_rank, contam_rank


def load_decontx_corrected(condition):
    out_dir = decontx_dir / condition
    input_dir = decontx_input_dir / condition

    mtx_path = out_dir / "decontx_corrected_counts_genes_by_cells.mtx"

    genes_path = out_dir / "genes.csv"
    if not genes_path.exists():
        genes_path = input_dir / "genes.csv"

    barcodes_path = out_dir / "barcodes.csv"
    if not barcodes_path.exists():
        barcodes_path = input_dir / "barcodes.csv"

    if not mtx_path.exists():
        raise FileNotFoundError(f"Missing DecontX corrected matrix: {mtx_path}")

    corrected_genes_by_cells = mmread(str(mtx_path)).tocsr()
    corrected_cells_by_genes = corrected_genes_by_cells.T.tocsr()

    genes = pd.read_csv(genes_path)["gene"].astype(str).values
    barcodes = pd.read_csv(barcodes_path)["barcode"].astype(str).values

    adata = sc.AnnData(
        X=corrected_cells_by_genes,
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes),
    )

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    return adata


def run_standard_umap_pipeline(adata, method_name):
    """
    Input should be raw count matrix.
    This runs the same downstream analysis for contaminated and corrected data.
    """
    adata = adata.copy()

    adata.layers["counts"] = as_csr(adata.X).copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat"
    )

    sc.pp.pca(
        adata,
        use_highly_variable=True,
        svd_solver="arpack"
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=10,
        n_pcs=40
    )

    sc.tl.umap(adata, random_state=12345)

    sc.tl.leiden(
        adata,
        resolution=0.5,
        key_added="leiden",
        random_state=12345
    )

    adata.obs["method"] = method_name

    return adata


def plot_three_umaps(condition, clean, contaminated, decontx, color_key, output_name, title_suffix):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    sc.pl.umap(
        clean,
        color=color_key,
        ax=axes[0],
        show=False,
        title=f"Clean\n{title_suffix}",
        frameon=False,
    )

    sc.pl.umap(
        contaminated,
        color=color_key,
        ax=axes[1],
        show=False,
        title=f"Contaminated\n{title_suffix}",
        frameon=False,
    )

    sc.pl.umap(
        decontx,
        color=color_key,
        ax=axes[2],
        show=False,
        title=f"DecontX corrected\n{title_suffix}",
        frameon=False,
    )

    plt.suptitle(f"UMAP comparison: {condition}", y=1.05, fontsize=14)
    plt.tight_layout()

    out_path = figure_dir / output_name
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved figure:", out_path)


# -----------------------------
# Main
# -----------------------------
print("Loading clean baseline...")
adata_clean = sc.read_h5ad(clean_path)

if "leiden" not in adata_clean.obs.columns:
    raise ValueError("Clean baseline does not contain Leiden labels. Please rerun 01_qc_preprocess.py.")

if "X_umap" not in adata_clean.obsm:
    raise ValueError("Clean baseline does not contain UMAP coordinates. Please rerun 01_qc_preprocess.py.")

# Default: only use two representative high-contamination conditions
if len(sys.argv) > 1:
    conditions = sys.argv[1:]
else:
    conditions = ["1gene_30pct", "100genes_30pct"]

print("Conditions to analyze:")
for c in conditions:
    print(" -", c)

records = []

for condition in conditions:
    print("\n==============================")
    print("Processing condition:", condition)
    print("==============================")

    contam_path = sim_dir / f"pbmc10k_contaminated_{condition}.h5ad"

    if not contam_path.exists():
        print("Skipping. Missing contaminated file:")
        print(contam_path)
        continue

    adata_contam = sc.read_h5ad(contam_path)
    adata_decontx = load_decontx_corrected(condition)

    common_cells = (
        adata_clean.obs_names
        .intersection(adata_contam.obs_names)
        .intersection(adata_decontx.obs_names)
    )

    common_genes = (
        adata_clean.var_names
        .intersection(adata_contam.var_names)
        .intersection(adata_decontx.var_names)
    )

    print("Common cells:", len(common_cells))
    print("Common genes:", len(common_genes))

    clean_sub = adata_clean[common_cells, common_genes].copy()
    contam_sub = adata_contam[common_cells, common_genes].copy()
    decontx_sub = adata_decontx[common_cells, common_genes].copy()

    # Use raw counts for contaminated and corrected data
    contam_counts = get_counts(contam_sub)
    decontx_counts = as_csr(decontx_sub.X)

    contam_for_umap = sc.AnnData(
        X=contam_counts,
        obs=clean_sub.obs.copy(),
        var=clean_sub.var.copy()
    )

    decontx_for_umap = sc.AnnData(
        X=decontx_counts,
        obs=clean_sub.obs.copy(),
        var=clean_sub.var.copy()
    )

    # Run identical UMAP/clustering pipeline
    contam_for_umap = run_standard_umap_pipeline(contam_for_umap, "No correction")
    decontx_for_umap = run_standard_umap_pipeline(decontx_for_umap, "DecontX")

    # Add clean Leiden labels to every object for consistent biological comparison
    clean_labels = clean_sub.obs["leiden"].astype("category")
    categories = clean_labels.cat.categories

    clean_sub.obs["clean_leiden"] = pd.Categorical(clean_labels.values, categories=categories)
    contam_for_umap.obs["clean_leiden"] = pd.Categorical(clean_labels.values, categories=categories)
    decontx_for_umap.obs["clean_leiden"] = pd.Categorical(clean_labels.values, categories=categories)

    # Quantitative clustering comparison
    clean_cluster = clean_sub.obs["leiden"].astype(str).values
    contam_cluster = contam_for_umap.obs["leiden"].astype(str).values
    decontx_cluster = decontx_for_umap.obs["leiden"].astype(str).values

    for method, cluster in [
        ("No correction", contam_cluster),
        ("DecontX", decontx_cluster),
    ]:
        ari = adjusted_rand_score(clean_cluster, cluster)
        nmi = normalized_mutual_info_score(clean_cluster, cluster)

        records.append({
            "condition": condition,
            "method": method,
            "n_cells": len(common_cells),
            "n_genes": len(common_genes),
            "clean_n_clusters": len(np.unique(clean_cluster)),
            "method_n_clusters": len(np.unique(cluster)),
            "ARI_vs_clean_leiden": ari,
            "NMI_vs_clean_leiden": nmi,
        })

    # Figure 1: color by original clean cluster labels
    plot_three_umaps(
        condition=condition,
        clean=clean_sub,
        contaminated=contam_for_umap,
        decontx=decontx_for_umap,
        color_key="clean_leiden",
        output_name=f"pbmc10k_{condition}_umap_by_clean_cluster.png",
        title_suffix="colored by clean cluster"
    )

    # Figure 2: color by each dataset's own Leiden clustering
    plot_three_umaps(
        condition=condition,
        clean=clean_sub,
        contaminated=contam_for_umap,
        decontx=decontx_for_umap,
        color_key="leiden",
        output_name=f"pbmc10k_{condition}_umap_by_own_leiden.png",
        title_suffix="colored by own Leiden cluster"
    )


summary = pd.DataFrame(records)

if len(summary) > 0:
    summary = summary.sort_values(
        by=["condition", "method"],
        key=lambda col: col.map(condition_order_key) if col.name == "condition" else col
    )

out_csv = table_dir / "pbmc10k_decontx_umap_clustering_impact.csv"
summary.to_csv(out_csv, index=False)

print("\nSaved clustering impact table:")
print(out_csv)

print("\nSummary:")
print(summary)

print("\nDone.")
