from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# -----------------------------
# 1. Set paths
# -----------------------------
base_dir = Path(__file__).resolve().parents[1]

input_path = base_dir / "data" / "processed" / "pbmc10k_clean_baseline.h5ad"
simulated_dir = base_dir / "data" / "simulated"
simulated_dir.mkdir(parents=True, exist_ok=True)

# -----------------------------
# 2. Load clean baseline
# -----------------------------
print("Loading clean baseline...")
adata = sc.read_h5ad(input_path)

if "counts" not in adata.layers:
    raise ValueError("No 'counts' layer found. Please check 01_qc_preprocess.py.")

counts = adata.layers["counts"]

if not sparse.issparse(counts):
    counts = sparse.csr_matrix(counts)
else:
    counts = counts.tocsr()

print("Clean baseline:")
print(adata)
print("Counts matrix shape:", counts.shape)

# -----------------------------
# 3. Choose ambient RNA genes
# -----------------------------
# We choose highly expressed non-mitochondrial, non-ribosomal genes
# to make the simulated ambient RNA signal visible and reproducible.

gene_names = np.array(adata.var_names)

total_gene_counts = np.asarray(counts.sum(axis=0)).ravel()

exclude = (
    pd.Series(gene_names).str.upper().str.startswith("MT-").values
    | pd.Series(gene_names).str.upper().str.startswith("RPS").values
    | pd.Series(gene_names).str.upper().str.startswith("RPL").values
)

valid_gene_idx = np.where((~exclude) & (total_gene_counts > 0))[0]

# Sort valid genes by total expression, highest first
valid_gene_idx_sorted = valid_gene_idx[np.argsort(total_gene_counts[valid_gene_idx])[::-1]]

ambient_1gene_idx = valid_gene_idx_sorted[:1]
ambient_100genes_idx = valid_gene_idx_sorted[:100]

print("\nSelected 1-gene ambient RNA:")
print(gene_names[ambient_1gene_idx])

print("\nSelected first 20 of 100 ambient RNA genes:")
print(gene_names[ambient_100genes_idx[:20]])

# Save ambient gene list
ambient_gene_df = pd.DataFrame({
    "ambient_set": ["1_gene"] * len(ambient_1gene_idx) + ["100_genes"] * len(ambient_100genes_idx),
    "gene": list(gene_names[ambient_1gene_idx]) + list(gene_names[ambient_100genes_idx]),
    "gene_index": list(ambient_1gene_idx) + list(ambient_100genes_idx),
    "total_clean_counts": list(total_gene_counts[ambient_1gene_idx]) + list(total_gene_counts[ambient_100genes_idx]),
})

ambient_gene_df.to_csv(simulated_dir / "simulated_ambient_genes.csv", index=False)

# -----------------------------
# 4. Function to simulate ambient RNA
# -----------------------------
def simulate_ambient_counts(clean_counts, ambient_gene_idx, contamination_fraction, random_seed=0):
    """
    Simulate ambient RNA contamination.

    contamination_fraction means:
        ambient_counts / (clean_counts + ambient_counts)

    Example:
        contamination_fraction = 0.15 means 15% of total UMIs after contamination
        come from simulated ambient RNA.
    """
    rng = np.random.default_rng(random_seed)

    n_cells, n_genes = clean_counts.shape

    clean_total_per_cell = np.asarray(clean_counts.sum(axis=1)).ravel()

    # If contamination_fraction = f,
    # ambient_total = clean_total * f / (1 - f)
    ambient_total_per_cell = np.rint(
        clean_total_per_cell * contamination_fraction / (1 - contamination_fraction)
    ).astype(int)

    # Ambient profile based on total expression of selected genes
    selected_gene_totals = np.asarray(clean_counts[:, ambient_gene_idx].sum(axis=0)).ravel()
    selected_gene_totals = selected_gene_totals.astype(float)

    if selected_gene_totals.sum() == 0:
        raise ValueError("Selected ambient genes have zero total counts.")

    ambient_prob = selected_gene_totals / selected_gene_totals.sum()

    rows = []
    cols = []
    vals = []

    for cell_idx, n_ambient in enumerate(ambient_total_per_cell):
        if n_ambient <= 0:
            continue

        added = rng.multinomial(n_ambient, ambient_prob)

        nonzero_local = np.where(added > 0)[0]
        if len(nonzero_local) == 0:
            continue

        rows.extend([cell_idx] * len(nonzero_local))
        cols.extend(ambient_gene_idx[nonzero_local])
        vals.extend(added[nonzero_local])

    ambient_matrix = sparse.csr_matrix(
        (vals, (rows, cols)),
        shape=clean_counts.shape,
        dtype=clean_counts.dtype,
    )

    contaminated_counts = clean_counts + ambient_matrix

    return contaminated_counts.tocsr(), ambient_matrix.tocsr(), ambient_total_per_cell


# -----------------------------
# 5. Simulate six conditions
# -----------------------------
conditions = []

ambient_sets = {
    "1gene": ambient_1gene_idx,
    "100genes": ambient_100genes_idx,
}

contamination_levels = [0.05, 0.15, 0.30]

summary_rows = []

for ambient_set_name, ambient_gene_idx in ambient_sets.items():
    for frac in contamination_levels:
        pct_label = int(frac * 100)
        condition_name = f"{ambient_set_name}_{pct_label}pct"

        print(f"\nSimulating condition: {condition_name}")

        contaminated_counts, ambient_added, ambient_total_per_cell = simulate_ambient_counts(
            clean_counts=counts,
            ambient_gene_idx=ambient_gene_idx,
            contamination_fraction=frac,
            random_seed=2026 + pct_label + len(ambient_gene_idx),
        )

        # Create new AnnData object for contaminated dataset
        adata_contam = sc.AnnData(
            X=contaminated_counts,
            obs=adata.obs.copy(),
            var=adata.var.copy(),
        )

        # Keep clean counts and ambient counts for ground-truth benchmarking
        adata_contam.layers["counts"] = contaminated_counts.copy()
        adata_contam.layers["clean_counts"] = counts.copy()
        adata_contam.layers["ambient_added"] = ambient_added.copy()

        # Keep clean clustering labels if available
        if "leiden" in adata.obs.columns:
            adata_contam.obs["clean_leiden"] = adata.obs["leiden"].astype(str).values

        adata_contam.obs["simulated_ambient_total"] = ambient_total_per_cell
        adata_contam.obs["target_contamination_fraction"] = frac

        adata_contam.uns["simulation"] = {
            "ambient_set": ambient_set_name,
            "n_ambient_genes": int(len(ambient_gene_idx)),
            "target_contamination_fraction": float(frac),
            "ambient_genes": list(gene_names[ambient_gene_idx]),
        }

        output_path = simulated_dir / f"pbmc10k_contaminated_{condition_name}.h5ad"
        adata_contam.write_h5ad(output_path)

        clean_total = counts.sum()
        ambient_total = ambient_added.sum()
        contaminated_total = contaminated_counts.sum()
        observed_frac = ambient_total / contaminated_total

        summary_rows.append({
            "condition": condition_name,
            "ambient_set": ambient_set_name,
            "n_ambient_genes": len(ambient_gene_idx),
            "target_contamination_fraction": frac,
            "total_clean_counts": int(clean_total),
            "total_ambient_added": int(ambient_total),
            "total_contaminated_counts": int(contaminated_total),
            "observed_contamination_fraction": float(observed_frac),
            "output_file": str(output_path.name),
        })

        print(f"Saved: {output_path.name}")
        print(f"Observed contamination fraction: {observed_frac:.4f}")

# -----------------------------
# 6. Save simulation summary
# -----------------------------
summary_df = pd.DataFrame(summary_rows)
summary_path = simulated_dir / "simulation_summary.csv"
summary_df.to_csv(summary_path, index=False)

print("\nSimulation summary:")
print(summary_df)

print("\nSaved summary to:")
print(summary_path)

print("\nDone.")
