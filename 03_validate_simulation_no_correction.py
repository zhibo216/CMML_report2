from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


# -----------------------------
# 1. Paths
# -----------------------------
base_dir = Path(__file__).resolve().parents[1]

clean_path = base_dir / "data" / "processed" / "pbmc10k_clean_baseline.h5ad"
sim_dir = base_dir / "data" / "simulated"
figure_dir = base_dir / "results" / "figures"
table_dir = base_dir / "results" / "tables"

figure_dir.mkdir(parents=True, exist_ok=True)
table_dir.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = figure_dir
sc.settings.verbosity = 2


# -----------------------------
# 2. Helper functions
# -----------------------------
def get_counts(adata):
    """
    Return raw count matrix.
    Prefer adata.layers['counts'] if available.
    """
    if "counts" in adata.layers:
        return adata.layers["counts"]
    return adata.X


def to_1d(x):
    return np.asarray(x).ravel()


def condition_from_filename(path):
    name = path.stem
    return name.replace("pbmc10k_contaminated_", "")


def plot_bar(df, x, y, title, ylabel, output_name):
    plt.figure(figsize=(8, 4))
    plt.bar(df[x], df[y])
    plt.xticks(rotation=45, ha="right")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(figure_dir / output_name, dpi=300)
    plt.close()


# -----------------------------
# 3. Load clean baseline
# -----------------------------
print("Loading clean baseline...")
adata_clean = sc.read_h5ad(clean_path)

if "leiden" not in adata_clean.obs:
    raise ValueError("Clean baseline does not contain adata.obs['leiden'].")

clean_counts = get_counts(adata_clean)

print("Clean baseline:")
print(adata_clean)


# -----------------------------
# 4. Find simulated files
# -----------------------------
sim_files = sorted(sim_dir.glob("pbmc10k_contaminated_*.h5ad"))

if len(sim_files) == 0:
    raise FileNotFoundError(f"No simulated h5ad files found in: {sim_dir}")

print("\nFound simulated files:")
for f in sim_files:
    print(f.name)


# -----------------------------
# 5. Benchmark contaminated data before correction
# -----------------------------
records = []

for sim_file in sim_files:
    condition = condition_from_filename(sim_file)
    print(f"\nProcessing condition: {condition}")

    adata_contam = sc.read_h5ad(sim_file)

    # Align cells and genes
    common_cells = adata_clean.obs_names.intersection(adata_contam.obs_names)
    common_genes = adata_clean.var_names.intersection(adata_contam.var_names)

    clean_sub = adata_clean[common_cells, common_genes].copy()
    contam_sub = adata_contam[common_cells, common_genes].copy()

    clean_mat = get_counts(clean_sub)
    contam_mat = get_counts(contam_sub)

    # Library size based contamination estimate
    clean_total_per_cell = to_1d(clean_mat.sum(axis=1))
    contam_total_per_cell = to_1d(contam_mat.sum(axis=1))
    added_total_per_cell = contam_total_per_cell - clean_total_per_cell

    observed_fraction = added_total_per_cell.sum() / contam_total_per_cell.sum()

    # Pseudobulk gene expression
    clean_gene_sum = to_1d(clean_mat.sum(axis=0))
    contam_gene_sum = to_1d(contam_mat.sum(axis=0))

    pearson = np.corrcoef(clean_gene_sum, contam_gene_sum)[0, 1]
    spearman = spearmanr(clean_gene_sum, contam_gene_sum).correlation

    clean_gene_mean = clean_gene_sum / clean_sub.n_obs
    contam_gene_mean = contam_gene_sum / contam_sub.n_obs

    mae_gene_mean = np.mean(np.abs(contam_gene_mean - clean_gene_mean))
    mse_gene_mean = np.mean((contam_gene_mean - clean_gene_mean) ** 2)

    # Run clustering on contaminated data
    adata_tmp = contam_sub.copy()
    adata_tmp.X = get_counts(adata_tmp).copy()
    adata_tmp.layers["counts"] = adata_tmp.X.copy()

    sc.pp.normalize_total(adata_tmp, target_sum=1e4)
    sc.pp.log1p(adata_tmp)
    sc.pp.highly_variable_genes(adata_tmp, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(adata_tmp, use_highly_variable=True, svd_solver="arpack")
    sc.pp.neighbors(adata_tmp, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata_tmp)
    sc.tl.leiden(adata_tmp, resolution=0.5, key_added="leiden")

    clean_labels = clean_sub.obs["leiden"].astype(str).values
    contam_labels = adata_tmp.obs["leiden"].astype(str).values

    ari = adjusted_rand_score(clean_labels, contam_labels)
    nmi = normalized_mutual_info_score(clean_labels, contam_labels)

    # Save UMAP for each contaminated condition
    sc.pl.umap(
        adata_tmp,
        color=["leiden"],
        title=f"Contaminated: {condition}",
        save=f"_contaminated_{condition}.png",
        show=False,
    )

    records.append(
        {
            "condition": condition,
            "n_cells": clean_sub.n_obs,
            "n_genes": clean_sub.n_vars,
            "observed_contamination_fraction": observed_fraction,
            "pseudobulk_pearson_clean_vs_contaminated": pearson,
            "pseudobulk_spearman_clean_vs_contaminated": spearman,
            "gene_mean_MAE_clean_vs_contaminated": mae_gene_mean,
            "gene_mean_MSE_clean_vs_contaminated": mse_gene_mean,
            "ARI_clean_leiden_vs_contaminated_leiden": ari,
            "NMI_clean_leiden_vs_contaminated_leiden": nmi,
        }
    )


# -----------------------------
# 6. Save summary table
# -----------------------------
summary_df = pd.DataFrame(records)
summary_path = table_dir / "pbmc10k_no_correction_simulation_benchmark.csv"
summary_df.to_csv(summary_path, index=False)

print("\nNo-correction benchmark summary:")
print(summary_df)

print("\nSaved table to:")
print(summary_path)


# -----------------------------
# 7. Save benchmark plots
# -----------------------------
plot_bar(
    summary_df,
    x="condition",
    y="observed_contamination_fraction",
    title="Observed contamination fraction",
    ylabel="Observed contamination fraction",
    output_name="pbmc10k_observed_contamination_fraction.png",
)

plot_bar(
    summary_df,
    x="condition",
    y="pseudobulk_pearson_clean_vs_contaminated",
    title="Pseudobulk correlation: clean vs contaminated",
    ylabel="Pearson correlation",
    output_name="pbmc10k_clean_vs_contaminated_pearson.png",
)

plot_bar(
    summary_df,
    x="condition",
    y="gene_mean_MSE_clean_vs_contaminated",
    title="Gene-level MSE: clean vs contaminated",
    ylabel="MSE",
    output_name="pbmc10k_clean_vs_contaminated_mse.png",
)

plot_bar(
    summary_df,
    x="condition",
    y="ARI_clean_leiden_vs_contaminated_leiden",
    title="Clustering similarity: clean vs contaminated",
    ylabel="Adjusted Rand Index",
    output_name="pbmc10k_clean_vs_contaminated_ari.png",
)

plot_bar(
    summary_df,
    x="condition",
    y="NMI_clean_leiden_vs_contaminated_leiden",
    title="Clustering similarity: clean vs contaminated",
    ylabel="Normalized Mutual Information",
    output_name="pbmc10k_clean_vs_contaminated_nmi.png",
)

print("\nSaved figures to:")
print(figure_dir)

print("\nDone.")
