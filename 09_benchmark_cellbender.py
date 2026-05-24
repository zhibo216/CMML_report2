from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Paths
# -----------------------------
base_dir = Path(__file__).resolve().parents[1]

clean_path = base_dir / "data" / "processed" / "pbmc10k_clean_baseline.h5ad"

cellbender_path = (
    base_dir
    / "results"
    / "cellbender_final"
    / "pbmc10k_cellbender_filtered.h5"
)

metrics_path = (
    base_dir
    / "results"
    / "cellbender_final"
    / "pbmc10k_cellbender_metrics.csv"
)

table_dir = base_dir / "results" / "tables"
figure_dir = base_dir / "results" / "figures"

table_dir.mkdir(parents=True, exist_ok=True)
figure_dir.mkdir(parents=True, exist_ok=True)

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


def expression_metrics(reference_mat, test_mat):
    reference_gene_sum = np.asarray(reference_mat.sum(axis=0)).ravel()
    test_gene_sum = np.asarray(test_mat.sum(axis=0)).ravel()

    pearson = np.corrcoef(reference_gene_sum, test_gene_sum)[0, 1]
    spearman = spearmanr(reference_gene_sum, test_gene_sum).correlation

    reference_gene_mean = reference_gene_sum / reference_mat.shape[0]
    test_gene_mean = test_gene_sum / test_mat.shape[0]

    mae = np.mean(np.abs(test_gene_mean - reference_gene_mean))
    mse = np.mean((test_gene_mean - reference_gene_mean) ** 2)

    return pearson, spearman, mae, mse


def run_scanpy_pipeline(adata, method_name):
    adata = adata.copy()
    adata.layers["counts"] = as_csr(adata.X).copy()

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=2000,
        flavor="seurat",
    )

    sc.pp.pca(
        adata,
        use_highly_variable=True,
        svd_solver="arpack",
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=10,
        n_pcs=40,
    )

    sc.tl.umap(adata, random_state=12345)

    sc.tl.leiden(
        adata,
        resolution=0.5,
        key_added="leiden",
        random_state=12345,
    )

    adata.obs["method"] = method_name

    return adata


def plot_umap_pair(adata_a, adata_b, color_key, title_a, title_b, out_name):
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))

    sc.pl.umap(
        adata_a,
        color=color_key,
        ax=axes[0],
        show=False,
        title=title_a,
        frameon=False,
    )

    sc.pl.umap(
        adata_b,
        color=color_key,
        ax=axes[1],
        show=False,
        title=title_b,
        frameon=False,
    )

    plt.tight_layout()

    out_path = figure_dir / out_name
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved figure:", out_path)


# -----------------------------
# Main
# -----------------------------
print("Checking input files...")
print("Clean baseline:", clean_path)
print("CellBender output:", cellbender_path)

if not clean_path.exists():
    raise FileNotFoundError(f"Missing clean baseline: {clean_path}")

if not cellbender_path.exists():
    raise FileNotFoundError(f"Missing CellBender output: {cellbender_path}")

print("\nLoading clean Cell Ranger filtered baseline...")
adata_clean = sc.read_h5ad(clean_path)
adata_clean.var_names_make_unique()
adata_clean.obs_names_make_unique()

print(adata_clean)

print("\nLoading CellBender corrected filtered matrix...")
adata_cb = sc.read_10x_h5(cellbender_path)
adata_cb.var_names_make_unique()
adata_cb.obs_names_make_unique()

print(adata_cb)

# Align cells and genes
common_cells = adata_clean.obs_names.intersection(adata_cb.obs_names)
common_genes = adata_clean.var_names.intersection(adata_cb.var_names)

print("\nCommon cells:", len(common_cells))
print("Common genes:", len(common_genes))

if len(common_cells) == 0:
    raise ValueError("No common cell barcodes found between Cell Ranger and CellBender outputs.")

if len(common_genes) == 0:
    raise ValueError("No common genes found between Cell Ranger and CellBender outputs.")

clean_sub = adata_clean[common_cells, common_genes].copy()
cb_sub = adata_cb[common_cells, common_genes].copy()

clean_counts = get_counts(clean_sub)
cb_counts = as_csr(cb_sub.X)

# -----------------------------
# Expression-level benchmark
# -----------------------------
pearson, spearman, mae, mse = expression_metrics(clean_counts, cb_counts)

clean_umi_per_cell = np.asarray(clean_counts.sum(axis=1)).ravel()
cb_umi_per_cell = np.asarray(cb_counts.sum(axis=1)).ravel()

clean_genes_per_cell = np.asarray((clean_counts > 0).sum(axis=1)).ravel()
cb_genes_per_cell = np.asarray((cb_counts > 0).sum(axis=1)).ravel()

summary = pd.DataFrame([
    {
        "comparison": "CellBender corrected vs Cell Ranger filtered baseline",
        "n_cellranger_cells": adata_clean.n_obs,
        "n_cellbender_cells": adata_cb.n_obs,
        "n_common_cells": len(common_cells),
        "n_common_genes": len(common_genes),
        "cell_overlap_fraction_vs_cellranger": len(common_cells) / adata_clean.n_obs,
        "cell_overlap_fraction_vs_cellbender": len(common_cells) / adata_cb.n_obs,
        "cellranger_total_umis_common_cells": clean_umi_per_cell.sum(),
        "cellbender_total_umis_common_cells": cb_umi_per_cell.sum(),
        "fraction_umis_remaining_after_cellbender": cb_umi_per_cell.sum() / clean_umi_per_cell.sum(),
        "median_umis_per_cell_cellranger": np.median(clean_umi_per_cell),
        "median_umis_per_cell_cellbender": np.median(cb_umi_per_cell),
        "median_genes_per_cell_cellranger": np.median(clean_genes_per_cell),
        "median_genes_per_cell_cellbender": np.median(cb_genes_per_cell),
        "pseudobulk_pearson_vs_cellranger": pearson,
        "pseudobulk_spearman_vs_cellranger": spearman,
        "gene_mean_MAE_vs_cellranger": mae,
        "gene_mean_MSE_vs_cellranger": mse,
    }
])

summary_path = table_dir / "pbmc10k_cellbender_benchmark.csv"
summary.to_csv(summary_path, index=False)

print("\nSaved CellBender benchmark table:")
print(summary_path)
print(summary.T)

# -----------------------------
# UMI / gene count comparison figure
# -----------------------------
qc_df = pd.DataFrame({
    "Cell Ranger filtered": clean_umi_per_cell,
    "CellBender corrected": cb_umi_per_cell,
})

plt.figure(figsize=(6, 4.5))
plt.boxplot(
    [clean_umi_per_cell, cb_umi_per_cell],
    labels=["Cell Ranger\nfiltered", "CellBender\ncorrected"],
    showfliers=False,
)
plt.ylabel("UMI counts per cell")
plt.title("UMI counts before and after CellBender")
plt.tight_layout()

out_path = figure_dir / "pbmc10k_cellbender_umi_count_comparison.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close()
print("Saved figure:", out_path)

plt.figure(figsize=(6, 4.5))
plt.boxplot(
    [clean_genes_per_cell, cb_genes_per_cell],
    labels=["Cell Ranger\nfiltered", "CellBender\ncorrected"],
    showfliers=False,
)
plt.ylabel("Detected genes per cell")
plt.title("Detected genes before and after CellBender")
plt.tight_layout()

out_path = figure_dir / "pbmc10k_cellbender_gene_count_comparison.png"
plt.savefig(out_path, dpi=300, bbox_inches="tight")
plt.close()
print("Saved figure:", out_path)

# -----------------------------
# Run Scanpy pipelines
# -----------------------------
print("\nRunning Scanpy pipeline on Cell Ranger baseline common cells...")
adata_cr_proc = sc.AnnData(
    X=clean_counts,
    obs=clean_sub.obs.copy(),
    var=clean_sub.var.copy(),
)

adata_cr_proc = run_scanpy_pipeline(
    adata_cr_proc,
    method_name="Cell Ranger filtered",
)

print("\nRunning Scanpy pipeline on CellBender corrected common cells...")
adata_cb_proc = sc.AnnData(
    X=cb_counts,
    obs=clean_sub.obs.copy(),
    var=clean_sub.var.copy(),
)

adata_cb_proc = run_scanpy_pipeline(
    adata_cb_proc,
    method_name="CellBender corrected",
)

# Add Cell Ranger cluster labels to CellBender object
cr_labels = adata_cr_proc.obs["leiden"].astype("category")
categories = cr_labels.cat.categories

adata_cr_proc.obs["cellranger_leiden"] = pd.Categorical(
    cr_labels.values,
    categories=categories,
)

adata_cb_proc.obs["cellranger_leiden"] = pd.Categorical(
    cr_labels.values,
    categories=categories,
)

# Cluster comparison
cr_cluster = adata_cr_proc.obs["leiden"].astype(str).values
cb_cluster = adata_cb_proc.obs["leiden"].astype(str).values

ari = adjusted_rand_score(cr_cluster, cb_cluster)
nmi = normalized_mutual_info_score(cr_cluster, cb_cluster)

cluster_summary = pd.DataFrame([
    {
        "comparison": "CellBender Leiden vs Cell Ranger Leiden",
        "n_common_cells": len(common_cells),
        "cellranger_n_clusters": len(np.unique(cr_cluster)),
        "cellbender_n_clusters": len(np.unique(cb_cluster)),
        "ARI_vs_cellranger_leiden": ari,
        "NMI_vs_cellranger_leiden": nmi,
    }
])

cluster_path = table_dir / "pbmc10k_cellbender_clustering_impact.csv"
cluster_summary.to_csv(cluster_path, index=False)

print("\nSaved clustering impact table:")
print(cluster_path)
print(cluster_summary.T)

# UMAP figures
plot_umap_pair(
    adata_cr_proc,
    adata_cb_proc,
    color_key="leiden",
    title_a="Cell Ranger filtered\nown Leiden",
    title_b="CellBender corrected\nown Leiden",
    out_name="pbmc10k_cellranger_vs_cellbender_umap_by_own_leiden.png",
)

plot_umap_pair(
    adata_cr_proc,
    adata_cb_proc,
    color_key="cellranger_leiden",
    title_a="Cell Ranger filtered\nCell Ranger clusters",
    title_b="CellBender corrected\ncolored by Cell Ranger clusters",
    out_name="pbmc10k_cellranger_vs_cellbender_umap_by_cellranger_cluster.png",
)

# -----------------------------
# Marker gene dotplots
# -----------------------------
marker_genes = [
    "CD3D", "CD3E", "IL7R",
    "MS4A1", "CD79A",
    "LYZ", "S100A8", "S100A9",
    "NKG7", "GNLY",
    "FCER1A", "CST3",
    "PPBP", "PF4",
]

marker_genes_present = [g for g in marker_genes if g in adata_cr_proc.var_names]

if len(marker_genes_present) > 0:
    sc.pl.dotplot(
        adata_cr_proc,
        marker_genes_present,
        groupby="leiden",
        save="_pbmc10k_cellranger_marker_dotplot.png",
        show=False,
    )

    sc.pl.dotplot(
        adata_cb_proc,
        marker_genes_present,
        groupby="leiden",
        save="_pbmc10k_cellbender_marker_dotplot.png",
        show=False,
    )

    print("Saved marker dotplots.")
else:
    print("No marker genes found for dotplot.")

# -----------------------------
# Save processed objects
# -----------------------------
adata_cr_proc.write_h5ad(
    base_dir / "data" / "processed" / "pbmc10k_cellranger_common_processed.h5ad"
)

adata_cb_proc.write_h5ad(
    base_dir / "data" / "processed" / "pbmc10k_cellbender_processed.h5ad"
)

print("\nSaved processed AnnData objects.")
print("\nDone.")
