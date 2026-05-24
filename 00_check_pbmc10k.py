from pathlib import Path
import scanpy as sc

base_dir = Path(__file__).resolve().parents[1]

raw_path = base_dir / "data" / "raw" / "pbmc10k" / "10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5"
filtered_path = base_dir / "data" / "raw" / "pbmc10k" / "10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5"

print("Checking files...")
print("Raw file exists:", raw_path.exists())
print("Filtered file exists:", filtered_path.exists())

print("\nReading filtered matrix...")
adata_filtered = sc.read_10x_h5(filtered_path)
adata_filtered.var_names_make_unique()
print(adata_filtered)

print("\nReading raw matrix...")
adata_raw = sc.read_10x_h5(raw_path)
adata_raw.var_names_make_unique()
print(adata_raw)

print("\nDone. Both files can be read by Scanpy.")
