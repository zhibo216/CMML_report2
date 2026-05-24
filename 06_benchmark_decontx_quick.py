from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.io import mmread
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt



base_dir = Path(__file__).resolve().parents[1]

clean_path = base_dir / "data" / "processed" / "pbmc10k_clean_baseline.h5ad"
sim_dir = base_dir / "data" / "simulated"
decontx_dir = base_dir / "data" / "decontx_output"
decontx_input_dir = base_dir / "data" / "decontx_input"

table_dir = base_dir / "results" / "tables"
figure_dir = base_dir / "results" / "figures"

table_dir.mkdir(parents=True, exist_ok=True)
figure_dir.mkdir(parents=True, exist_ok=True)


def as_csr(x):
    if sparse.issparse(x):
        return x.tocsr()
    return sparse.csr_matrix(x)


def get_counts(adata):
    if "counts" in adata.layers:
        return as_csr(adata.layers["counts"])
    return as_csr(adata.X)


def parse_contamination_level(condition):
    condition = str(condition)

    if "30pct" in condition:
        return 0.30
    if "15pct" in condition:
        return 0.15
    if "5pct" in condition:
        return 0.05

    return np.nan


def parse_ambient_set(condition):
    if condition.startswith("1gene"):
        return "1gene"
    if condition.startswith("100genes"):
        return "100genes"
    return "unknown"

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

    return adata


def expression_metrics(clean_mat, test_mat):
    clean_gene_sum = np.asarray(clean_mat.sum(axis=0)).ravel()
    test_gene_sum = np.asarray(test_mat.sum(axis=0)).ravel()

    pearson = np.corrcoef(clean_gene_sum, test_gene_sum)[0, 1]
    spearman = spearmanr(clean_gene_sum, test_gene_sum).correlation

    clean_gene_mean = clean_gene_sum / clean_mat.shape[0]
    test_gene_mean = test_gene_sum / test_mat.shape[0]

    mae = np.mean(np.abs(test_gene_mean - clean_gene_mean))
    mse = np.mean((test_gene_mean - clean_gene_mean) ** 2)

    return pearson, spearman, mae, mse


print("Loading clean baseline...")
adata_clean = sc.read_h5ad(clean_path)

conditions = sorted(
    [
        p.name for p in decontx_dir.iterdir()
        if p.is_dir() and (p / "decontx_corrected_counts_genes_by_cells.mtx").exists()
    ],
    key=condition_order_key
)


print("Found DecontX output conditions:")
for c in conditions:
    print(" -", c)

records = []

for condition in conditions:
    print(f"\nBenchmarking: {condition}")

    contam_path = sim_dir / f"pbmc10k_contaminated_{condition}.h5ad"

    if not contam_path.exists():
        print(f"Skipping {condition}, missing contaminated file: {contam_path}")
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

    clean_sub = adata_clean[common_cells, common_genes].copy()
    contam_sub = adata_contam[common_cells, common_genes].copy()
    decontx_sub = adata_decontx[common_cells, common_genes].copy()

    clean_counts = get_counts(clean_sub)
    contam_counts = get_counts(contam_sub)
    decontx_counts = as_csr(decontx_sub.X)

    for method, mat in [
        ("No correction", contam_counts),
        ("DecontX", decontx_counts),
    ]:
        pearson, spearman, mae, mse = expression_metrics(clean_counts, mat)

        mean_estimated_contam = np.nan
        median_estimated_contam = np.nan

        meta_path = decontx_dir / condition / "decontx_cell_metadata.csv"
        if method == "DecontX" and meta_path.exists():
            meta = pd.read_csv(meta_path)
            contam_cols = [c for c in meta.columns if "contamination" in c.lower()]
            if len(contam_cols) > 0:
                vals = meta[contam_cols[0]]
                mean_estimated_contam = vals.mean()
                median_estimated_contam = vals.median()

        records.append({
            "condition": condition,
            "ambient_set": parse_ambient_set(condition),
            "designed_contamination_fraction": parse_contamination_level(condition),
            "method": method,
            "n_cells": clean_sub.n_obs,
            "n_genes": clean_sub.n_vars,
            "pseudobulk_pearson_vs_clean": pearson,
            "pseudobulk_spearman_vs_clean": spearman,
            "gene_mean_MAE_vs_clean": mae,
            "gene_mean_MSE_vs_clean": mse,
            "mean_decontx_estimated_contamination": mean_estimated_contam,
            "median_decontx_estimated_contamination": median_estimated_contam,
        })


summary = pd.DataFrame(records)

out_csv = table_dir / "pbmc10k_decontx_quick_benchmark.csv"
summary.to_csv(out_csv, index=False)

print("\nSaved benchmark table:")
print(out_csv)

print("\nSummary:")
print(summary)


def plot_metric(metric, ylabel, output_name):
    ordered_conditions = sorted(
        summary["condition"].unique(),
        key=condition_order_key
    )

    pivot = summary.pivot_table(
        index="condition",
        columns="method",
        values=metric,
        aggfunc="mean"
    )

    pivot = pivot.reindex(ordered_conditions)

    methods = [m for m in ["No correction", "DecontX"] if m in pivot.columns]

    ax = pivot[methods].plot(kind="bar", figsize=(10, 4.8))

    ax.set_ylabel(ylabel)
    ax.set_xlabel("Condition")
    ax.set_title(ylabel)

    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    out_path = figure_dir / output_name
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()

    print("Saved figure:", out_path)

plot_metric(
    "pseudobulk_pearson_vs_clean",
    "Pseudobulk Pearson correlation vs clean",
    "pbmc10k_decontx_quick_pearson.png",
)

plot_metric(
    "gene_mean_MSE_vs_clean",
    "Gene-level MSE vs clean",
    "pbmc10k_decontx_quick_mse.png",
)

plot_metric(
    "gene_mean_MAE_vs_clean",
    "Gene-level MAE vs clean",
    "pbmc10k_decontx_quick_mae.png",
)


print("\nSaved figures to:")
print(figure_dir)
print("\nDone.")
