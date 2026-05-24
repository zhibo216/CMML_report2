# PBMC10k Ambient RNA Contamination Benchmark

This repository contains the analysis workflow for **CMML3 Miniproject 8: Benchmarking ambient RNA detection and correction methods in single-cell RNA-seq data**. The project uses the 10x Genomics PBMC10k dataset to simulate ambient RNA contamination, evaluate DecontX on simulated called-cell matrices, and compare CellBender output with the Cell Ranger filtered baseline.

The scripts are intended to be run in numerical order. Earlier scripts generate intermediate files that are required by later steps.

## Project aims

This workflow addresses three main questions:

1. How do ambient RNA complexity and contamination level affect expression error and downstream clustering?
2. Can DecontX reduce simulated ambient RNA contamination when an operational clean reference is available?
3. Does CellBender preserve broad PBMC expression and clustering structure when applied to the raw droplet matrix?

DecontX and CellBender are evaluated in complementary settings rather than as a direct head-to-head comparison. DecontX is tested on simulated contaminated called-cell matrices with an operational clean reference. CellBender is applied to raw droplet data and compared with the Cell Ranger filtered baseline because it requires empty-droplet information.

## Recommended project structure

Place all scripts in a `scripts` directory under the project root. The scripts use the parent directory of `scripts` as the project root.

```text
project8_ambientRNA
├── scripts
│   ├── 00_check_pbmc10k.py
│   ├── 01_qc_preprocess.py
│   ├── 02_simulate_ambientRNA.py
│   ├── 03_validate_simulation_no_correction.py
│   ├── 04_export_decontx_inputs.py
│   ├── 05_run_decontx_optimized.R
│   ├── 06_benchmark_decontx_quick.py
│   ├── 07_plot_decontx_contamination_estimate.py
│   ├── 08_umap_clustering_impact.py
│   └── 09_benchmark_cellbender.py
├── data
│   ├── raw
│   │   └── pbmc10k
│   ├── processed
│   ├── simulated
│   ├── decontx_input
│   └── decontx_output
└── results
    ├── tables
    ├── figures
    └── cellbender_final
```

## Input data

Place the original PBMC10k files in:

```text
data/raw/pbmc10k
```

The directory should contain the following 10x HDF5 files:

```text
10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5
10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5
```

The filtered feature-barcode matrix is used to generate the quality-controlled operational clean reference, run downstream clustering, and simulate ambient RNA contamination. The raw feature-barcode matrix is required for the CellBender workflow.

If you want to run the CellBender comparison, prepare the following file before running script 09:

```text
results/cellbender_final/pbmc10k_cellbender_filtered.h5
```

## Software requirements

The Python workflow requires the following main packages:

```text
scanpy
numpy
pandas
scipy
scikit-learn
matplotlib
```

The R workflow for DecontX requires:

```text
Matrix
SingleCellExperiment
SummarizedExperiment
celda
```

DecontX is provided through the `celda` package. Make sure all R packages are installed before running the R script.

## Running the workflow

Run all commands from the project root unless stated otherwise.

### 1. Check PBMC10k input files

```bash
python scripts/00_check_pbmc10k.py
```

This script checks that the raw and filtered 10x HDF5 files exist and confirms that Scanpy can read them correctly.

### 2. Generate the operational clean reference

```bash
python scripts/01_qc_preprocess.py
```

This script reads the Cell Ranger filtered matrix, calculates quality-control metrics, filters low-quality cells and genes, and applies a standard Scanpy workflow. Raw counts are stored in `adata.layers["counts"]` before normalisation.

Main steps include:

1. Cell and gene filtering
2. Storage of raw counts
3. Total-count normalisation and log transformation
4. Selection of 2,000 highly variable genes
5. PCA, nearest-neighbour graph construction, UMAP and Leiden clustering
6. Marker-gene analysis

Main outputs:

```text
data/processed/pbmc10k_clean_baseline.h5ad
data/processed/pbmc10k_qc_metrics_before_filtering.csv
data/processed/pbmc10k_qc_metrics_after_filtering.csv
data/processed/pbmc10k_clean_baseline_marker_genes.csv
results/figures
```

The file `pbmc10k_clean_baseline.h5ad` is required for the downstream simulation and benchmarking steps.

### 3. Simulate ambient RNA contamination

```bash
python scripts/02_simulate_ambientRNA.py
```

This script simulates ambient RNA contamination using raw counts from the operational clean reference. The ambient profiles are based on highly expressed non-mitochondrial and non-ribosomal genes.

Six simulated conditions are generated:

```text
1gene_5pct
1gene_15pct
1gene_30pct
100genes_5pct
100genes_15pct
100genes_30pct
```

Each condition produces a contaminated `.h5ad` file and records the simulated ambient genes and contamination settings.

Main outputs:

```text
data/simulated/pbmc10k_contaminated_1gene_5pct.h5ad
data/simulated/pbmc10k_contaminated_1gene_15pct.h5ad
data/simulated/pbmc10k_contaminated_1gene_30pct.h5ad
data/simulated/pbmc10k_contaminated_100genes_5pct.h5ad
data/simulated/pbmc10k_contaminated_100genes_15pct.h5ad
data/simulated/pbmc10k_contaminated_100genes_30pct.h5ad
data/simulated/simulated_ambient_genes.csv
data/simulated/simulation_summary.csv
```

### 4. Evaluate the no-correction baseline

```bash
python scripts/03_validate_simulation_no_correction.py
```

This script compares each contaminated matrix with the operational clean reference. It is used to confirm that simulated ambient RNA changes expression and downstream structure before correction.

Metrics include:

1. Observed contamination fraction
2. Pseudobulk Pearson correlation
3. Pseudobulk Spearman correlation
4. Gene-level MAE and MSE
5. ARI and NMI between clean and contaminated Leiden clusters
6. UMAP visualisation for each simulated condition

Main outputs:

```text
results/tables/pbmc10k_no_correction_simulation_benchmark.csv
results/figures/pbmc10k_observed_contamination_fraction.png
results/figures/pbmc10k_clean_vs_contaminated_pearson.png
results/figures/pbmc10k_clean_vs_contaminated_mse.png
results/figures/pbmc10k_clean_vs_contaminated_ari.png
results/figures/pbmc10k_clean_vs_contaminated_nmi.png
```

### 5. Export DecontX input files

```bash
python scripts/04_export_decontx_inputs.py
```

DecontX is run in R, so this step converts each simulated `.h5ad` file into matrix and metadata files that are easier to load in R.

For each condition, the script exports:

```text
counts_genes_by_cells.mtx
genes.csv
barcodes.csv
clusters.csv
```

In this implementation, `clusters.csv` contains cluster labels generated from the contaminated matrix unless the script is modified to use alternative labels.

Main output directory:

```text
data/decontx_input
```

### 6. Run DecontX

```bash
Rscript scripts/05_run_decontx_optimized.R
```

By default, the script runs DecontX on all conditions in `data/decontx_input`.

To run only selected conditions, pass the condition names as command-line arguments:

```bash
Rscript scripts/05_run_decontx_optimized.R 1gene_30pct 100genes_30pct
```

Each condition produces corrected counts and cell-level DecontX contamination estimates.

Main output directory:

```text
data/decontx_output
```

Main output files for each condition:

```text
decontx_corrected_counts_genes_by_cells.mtx
decontx_cell_metadata.csv
genes.csv
barcodes.csv
```

### 7. Benchmark DecontX expression-level correction

```bash
python scripts/06_benchmark_decontx_quick.py
```

This script compares the no-correction and DecontX-corrected matrices against the operational clean reference.

Metrics include:

1. Pseudobulk Pearson correlation
2. Pseudobulk Spearman correlation
3. Gene-level MAE
4. Gene-level MSE
5. Mean and median DecontX-estimated contamination

Main outputs:

```text
results/tables/pbmc10k_decontx_quick_benchmark.csv
results/figures/pbmc10k_decontx_quick_pearson.png
results/figures/pbmc10k_decontx_quick_mse.png
results/figures/pbmc10k_decontx_quick_mae.png
```

### 8. Plot DecontX-estimated contamination

```bash
python scripts/07_plot_decontx_contamination_estimate.py
```

This script compares the designed contamination fractions with the mean and median contamination fractions estimated by DecontX.

Main output:

```text
results/figures/pbmc10k_decontx_estimated_vs_designed_contamination.png
```

### 9. Evaluate UMAP and clustering impact

```bash
python scripts/08_umap_clustering_impact.py
```

By default, this script analyses the two high-contamination conditions:

```text
1gene_30pct
100genes_30pct
```

Additional conditions can be passed manually:

```bash
python scripts/08_umap_clustering_impact.py 1gene_5pct 1gene_15pct 1gene_30pct 100genes_5pct 100genes_15pct 100genes_30pct
```

The script applies the same downstream Scanpy workflow to contaminated and DecontX-corrected matrices, then compares their Leiden clusters with the operational clean reference using ARI and NMI.

Main outputs:

```text
results/tables/pbmc10k_decontx_umap_clustering_impact.csv
results/figures/pbmc10k_1gene_30pct_umap_by_clean_cluster.png
results/figures/pbmc10k_1gene_30pct_umap_by_own_leiden.png
results/figures/pbmc10k_100genes_30pct_umap_by_clean_cluster.png
results/figures/pbmc10k_100genes_30pct_umap_by_own_leiden.png
```

If more conditions are provided, the script generates corresponding UMAP figures for each condition.

### 10. Run the CellBender comparison

```bash
python scripts/09_benchmark_cellbender.py
```

This optional step compares the CellBender-corrected matrix with the Cell Ranger filtered baseline. It requires the CellBender output file to be available before running the script:

```text
results/cellbender_final/pbmc10k_cellbender_filtered.h5
```

The comparison includes:

1. Common cells and common genes
2. UMI counts and detected genes per cell
3. Pseudobulk expression correlation
4. Gene-level MAE and MSE
5. Leiden clustering ARI and NMI
6. UMAP comparison
7. Marker-gene dot plot
8. Processed AnnData objects for downstream checks

Main outputs:

```text
results/tables/pbmc10k_cellbender_benchmark.csv
results/tables/pbmc10k_cellbender_clustering_impact.csv
results/figures/pbmc10k_cellbender_umi_count_comparison.png
results/figures/pbmc10k_cellbender_gene_count_comparison.png
results/figures/pbmc10k_cellranger_vs_cellbender_umap_by_own_leiden.png
results/figures/pbmc10k_cellranger_vs_cellbender_umap_by_cellranger_cluster.png
data/processed/pbmc10k_cellranger_common_processed.h5ad
data/processed/pbmc10k_cellbender_processed.h5ad
```

## Running the main pipeline

To run the DecontX simulation benchmark, use:

```bash
python scripts/00_check_pbmc10k.py
python scripts/01_qc_preprocess.py
python scripts/02_simulate_ambientRNA.py
python scripts/03_validate_simulation_no_correction.py
python scripts/04_export_decontx_inputs.py
Rscript scripts/05_run_decontx_optimized.R
python scripts/06_benchmark_decontx_quick.py
python scripts/07_plot_decontx_contamination_estimate.py
python scripts/08_umap_clustering_impact.py
```

Run the CellBender comparison separately after preparing the CellBender output:

```bash
python scripts/09_benchmark_cellbender.py
```

## Interpreting the outputs

Main result tables are saved in:

```text
results/tables
```

Main figures are saved in:

```text
results/figures
```

Key outputs to inspect first:

```text
results/tables/pbmc10k_no_correction_simulation_benchmark.csv
results/tables/pbmc10k_decontx_quick_benchmark.csv
results/tables/pbmc10k_decontx_umap_clustering_impact.csv
results/figures/pbmc10k_decontx_estimated_vs_designed_contamination.png
```

For DecontX, the main benchmark should be interpreted using multiple metrics. A decrease in gene-level MSE indicates reduced large gene-specific errors, but changes in MAE, pseudobulk correlation, ARI and NMI should also be inspected. The DecontX contamination estimate should be treated as a model-derived value rather than a calibrated recovery of the designed contamination fraction.

For CellBender, the comparison is made against the Cell Ranger filtered baseline rather than a true clean ground truth. The results should therefore be interpreted as evidence of background removal and structure preservation, not as a direct ranking against DecontX.

## Notes and troubleshooting

1. Keep the scripts inside the `scripts` directory. The Python scripts locate the project root by moving one level up from the script directory.
2. Run the R script from the project root so that it can find `data/decontx_input`.
3. Run `01_qc_preprocess.py` before any simulation step because later scripts require `data/processed/pbmc10k_clean_baseline.h5ad`.
4. Run `04_export_decontx_inputs.py` before `05_run_decontx_optimized.R`.
5. Run `06_benchmark_decontx_quick.py` before `07_plot_decontx_contamination_estimate.py`.
6. Run `09_benchmark_cellbender.py` only after the CellBender output file has been generated and placed in `results/cellbender_final`.
7. If a script cannot find an expected file, first check that all earlier scripts completed successfully and that the project directory structure matches the layout shown above.
