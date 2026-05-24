from pathlib import Path
import pandas as pd
import numpy as np

base_dir = Path(__file__).resolve().parents[1]
table_dir = base_dir / "results" / "tables"
figure_dir = base_dir / "results" / "figures"
cellbender_dir = base_dir / "results" / "cellbender_final"

out_dir = table_dir
out_dir.mkdir(parents=True, exist_ok=True)


def safe_read(path):
    if path.exists():
        return pd.read_csv(path)
    print(f"Missing file: {path}")
    return None


# -----------------------------
# 1. DecontX summary
# -----------------------------
decontx_path = table_dir / "pbmc10k_decontx_quick_benchmark.csv"
decontx = safe_read(decontx_path)

if decontx is not None:
    wide = decontx.pivot_table(
        index=["condition", "ambient_set", "designed_contamination_fraction"],
        columns="method",
        values=[
            "pseudobulk_pearson_vs_clean",
            "pseudobulk_spearman_vs_clean",
            "gene_mean_MAE_vs_clean",
            "gene_mean_MSE_vs_clean",
            "mean_decontx_estimated_contamination",
            "median_decontx_estimated_contamination",
        ],
        aggfunc="first"
    )

    wide.columns = ["__".join([str(x) for x in col if str(x) != ""]) for col in wide.columns]
    wide = wide.reset_index()

    if (
        "gene_mean_MSE_vs_clean__No correction" in wide.columns
        and "gene_mean_MSE_vs_clean__DecontX" in wide.columns
    ):
        wide["MSE_reduction_by_DecontX"] = (
            wide["gene_mean_MSE_vs_clean__No correction"]
            - wide["gene_mean_MSE_vs_clean__DecontX"]
        )

        wide["MSE_percent_reduction_by_DecontX"] = (
            wide["MSE_reduction_by_DecontX"]
            / wide["gene_mean_MSE_vs_clean__No correction"]
            * 100
        )

    if (
        "gene_mean_MAE_vs_clean__No correction" in wide.columns
        and "gene_mean_MAE_vs_clean__DecontX" in wide.columns
    ):
        wide["MAE_change_by_DecontX"] = (
            wide["gene_mean_MAE_vs_clean__DecontX"]
            - wide["gene_mean_MAE_vs_clean__No correction"]
        )

    out_path = out_dir / "report_decontx_summary.csv"
    wide.to_csv(out_path, index=False)
    print("Saved:", out_path)

    print("\nDecontX summary:")
    cols_to_show = [
        "condition",
        "gene_mean_MSE_vs_clean__No correction",
        "gene_mean_MSE_vs_clean__DecontX",
        "MSE_percent_reduction_by_DecontX",
        "gene_mean_MAE_vs_clean__No correction",
        "gene_mean_MAE_vs_clean__DecontX",
        "MAE_change_by_DecontX",
        "pseudobulk_pearson_vs_clean__No correction",
        "pseudobulk_pearson_vs_clean__DecontX",
    ]
    cols_to_show = [c for c in cols_to_show if c in wide.columns]
    print(wide[cols_to_show].to_string(index=False))


# -----------------------------
# 2. DecontX clustering summary
# -----------------------------
decontx_cluster_path = table_dir / "pbmc10k_decontx_umap_clustering_impact.csv"
decontx_cluster = safe_read(decontx_cluster_path)

if decontx_cluster is not None:
    out_path = out_dir / "report_decontx_clustering_summary.csv"
    decontx_cluster.to_csv(out_path, index=False)
    print("\nSaved:", out_path)
    print("\nDecontX clustering impact:")
    print(decontx_cluster.to_string(index=False))


# -----------------------------
# 3. CellBender summary
# -----------------------------
cellbender_path = table_dir / "pbmc10k_cellbender_benchmark.csv"
cellbender = safe_read(cellbender_path)

if cellbender is not None:
    out_path = out_dir / "report_cellbender_summary.csv"
    cellbender.to_csv(out_path, index=False)
    print("\nSaved:", out_path)

    print("\nCellBender benchmark summary:")
    print(cellbender.T.to_string())


cellbender_cluster_path = table_dir / "pbmc10k_cellbender_clustering_impact.csv"
cellbender_cluster = safe_read(cellbender_cluster_path)

if cellbender_cluster is not None:
    out_path = out_dir / "report_cellbender_clustering_summary.csv"
    cellbender_cluster.to_csv(out_path, index=False)
    print("\nSaved:", out_path)

    print("\nCellBender clustering impact:")
    print(cellbender_cluster.T.to_string())


# -----------------------------
# 4. CellBender metrics file
# -----------------------------
metrics_path = cellbender_dir / "pbmc10k_cellbender_metrics.csv"

if metrics_path.exists():
    metrics = pd.read_csv(metrics_path)
    out_path = out_dir / "report_cellbender_original_metrics.csv"
    metrics.to_csv(out_path, index=False)
    print("\nSaved:", out_path)
    print("\nCellBender original metrics columns:")
    print(list(metrics.columns))
else:
    print("\nCellBender metrics file not found:", metrics_path)


# -----------------------------
# 5. Figure checklist for report
# -----------------------------
figures = [
    {
        "section": "Data preprocessing",
        "figure": "umap_pbmc10k_clean_baseline.png",
        "use_in_report": "Show clean PBMC10k baseline UMAP and QC-related structure.",
        "priority": "High",
    },
    {
        "section": "Simulation effect",
        "figure": "pbmc10k_clean_vs_contaminated_mse.png",
        "use_in_report": "Show that simulated ambient RNA increases error relative to clean baseline.",
        "priority": "High",
    },
    {
        "section": "DecontX expression benchmark",
        "figure": "pbmc10k_decontx_quick_mse.png",
        "use_in_report": "Main DecontX result: correction reduces gene-level MSE.",
        "priority": "High",
    },
    {
        "section": "DecontX expression benchmark",
        "figure": "pbmc10k_decontx_quick_pearson.png",
        "use_in_report": "Show pseudobulk Pearson correlation changes little after correction.",
        "priority": "Medium",
    },
    {
        "section": "DecontX contamination estimate",
        "figure": "pbmc10k_decontx_estimated_vs_designed_contamination.png",
        "use_in_report": "Compare designed contamination fraction with DecontX-estimated contamination.",
        "priority": "High",
    },
    {
        "section": "DecontX UMAP impact",
        "figure": "pbmc10k_100genes_30pct_umap_by_clean_cluster.png",
        "use_in_report": "Show UMAP/clustering structure under high multi-gene contamination.",
        "priority": "High",
    },
    {
        "section": "CellBender count impact",
        "figure": "pbmc10k_cellbender_umi_count_comparison.png",
        "use_in_report": "Show UMI counts decrease after CellBender background removal.",
        "priority": "High",
    },
    {
        "section": "CellBender count impact",
        "figure": "pbmc10k_cellbender_gene_count_comparison.png",
        "use_in_report": "Show detected genes per cell before and after CellBender.",
        "priority": "Medium",
    },
    {
        "section": "CellBender UMAP impact",
        "figure": "pbmc10k_cellranger_vs_cellbender_umap_by_cellranger_cluster.png",
        "use_in_report": "Compare UMAP structure before and after CellBender using Cell Ranger cluster labels.",
        "priority": "High",
    },
    {
        "section": "CellBender marker genes",
        "figure": "dotplot_pbmc10k_cellbender_marker_dotplot.png",
        "use_in_report": "Show immune marker expression after CellBender correction.",
        "priority": "Medium",
    },
]

fig_df = pd.DataFrame(figures)
fig_df["exists"] = fig_df["figure"].apply(lambda x: (figure_dir / x).exists())

fig_path = out_dir / "report_figure_checklist.csv"
fig_df.to_csv(fig_path, index=False)

print("\nSaved:", fig_path)
print("\nFigure checklist:")
print(fig_df.to_string(index=False))


# -----------------------------
# 6. Plain text report notes
# -----------------------------
notes = []

notes.append("Project 8 report-ready summary")
notes.append("=" * 40)
notes.append("")
notes.append("Recommended report structure:")
notes.append("1. Introduction: ambient RNA and why correction matters.")
notes.append("2. Dataset: 10x PBMC10k raw and filtered matrices.")
notes.append("3. Preprocessing: QC, normalization, PCA, UMAP, Leiden clustering.")
notes.append("4. Simulation design: 1-gene vs 100-gene ambient profiles; 5%, 15%, 30% contamination.")
notes.append("5. Methods: No correction, DecontX, CellBender.")
notes.append("6. Metrics: Pearson/Spearman, MAE/MSE, ARI/NMI, UMI count changes, UMAPs, marker genes.")
notes.append("7. Results: DecontX simulation benchmark and CellBender raw-data benchmark.")
notes.append("8. Discussion: method strengths, limitations, and non-comparability of simulation vs raw-data workflows.")
notes.append("")
notes.append("Important interpretation point:")
notes.append("DecontX and CellBender are not evaluated under exactly the same input setting here.")
notes.append("DecontX is evaluated on simulated cell-level contaminated matrices with known clean ground truth.")
notes.append("CellBender is evaluated on the raw 10x droplet matrix, because it requires empty droplets/background barcodes.")
notes.append("Therefore, compare them conceptually and workflow-wise, not as a direct head-to-head numerical ranking.")
notes.append("")

notes_path = out_dir / "project8_report_notes.txt"
notes_path.write_text("\n".join(notes), encoding="utf-8")

print("\nSaved:", notes_path)
print("\nDone.")
