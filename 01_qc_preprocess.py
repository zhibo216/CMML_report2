from pathlib import Path
import scanpy as sc


base_dir = Path(__file__).resolve().parents[1]

filtered_path = (
    base_dir
    / "data"
    / "raw"
    / "pbmc10k"
    / "10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5"
)

processed_dir = base_dir / "data" / "processed"
figure_dir = base_dir / "results" / "figures"

processed_dir.mkdir(parents=True, exist_ok=True)
figure_dir.mkdir(parents=True, exist_ok=True)

sc.settings.figdir = figure_dir
sc.settings.verbosity = 3

print("Reading filtered PBMC10k matrix...")
adata = sc.read_10x_h5(filtered_path)
adata.var_names_make_unique()
adata.obs_names_make_unique()

print("Initial AnnData:")
print(adata)


adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

sc.pp.calculate_qc_metrics(
    adata,
    qc_vars=["mt"],
    percent_top=None,
    log1p=False,
    inplace=True,
)

adata.obs.to_csv(processed_dir / "pbmc10k_qc_metrics_before_filtering.csv")

print("\nQC summary before filtering:")
print(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_pbmc10k_qc_before_filtering.png",
    show=False,
)

print("\nFiltering cells...")

adata = adata[adata.obs["n_genes_by_counts"] >= 200, :].copy()
adata = adata[adata.obs["total_counts"] >= 500, :].copy()
adata = adata[adata.obs["pct_counts_mt"] < 20, :].copy()

sc.pp.filter_genes(adata, min_cells=3)

print("\nAfter filtering:")
print(adata)

print("\nQC summary after filtering:")
print(adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].describe())

adata.obs.to_csv(processed_dir / "pbmc10k_qc_metrics_after_filtering.csv")

# Save raw count matrix before normalization.
# This layer will be used later for ambient RNA simulation.
adata.layers["counts"] = adata.X.copy()


print("\nNormalizing and log-transforming...")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.raw = adata

print("\nSelecting highly variable genes...")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat",
)

print("\nRunning PCA...")
sc.pp.pca(
    adata,
    use_highly_variable=True,
    svd_solver="arpack",
)

print("\nRunning neighbors...")
sc.pp.neighbors(
    adata,
    n_neighbors=10,
    n_pcs=40,
)

print("\nRunning UMAP...")
sc.tl.umap(adata)

print("\nRunning Leiden clustering...")
sc.tl.leiden(
    adata,
    resolution=0.5,
    key_added="leiden",
)

sc.pl.umap(
    adata,
    color=["leiden", "total_counts", "n_genes_by_counts", "pct_counts_mt"],
    save="_pbmc10k_clean_baseline.png",
    show=False,
)

print("\nFinding marker genes...")

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
)

marker_df = sc.get.rank_genes_groups_df(adata, group=None)
marker_df.to_csv(
    processed_dir / "pbmc10k_clean_baseline_marker_genes.csv",
    index=False,
)

sc.pl.rank_genes_groups(
    adata,
    n_genes=5,
    sharey=False,
    save="_pbmc10k_clean_baseline_marker_genes.png",
    show=False,
)

output_path = processed_dir / "pbmc10k_clean_baseline.h5ad"
adata.write_h5ad(output_path)

print("\nSaved clean baseline to:")
print(output_path)

print("\nDone.")
