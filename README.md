# PBMC10k ambient RNA contamination benchmark

这个项目用 PBMC10k 单细胞数据做 ambient RNA 污染模拟，并比较未校正、DecontX 校正和 CellBender 结果之间的差异。代码整体按编号顺序运行，前面的脚本会生成后面脚本需要的中间文件。

项目主要关注三件事：

1. 在干净的 PBMC10k 数据上人为加入不同程度的 ambient RNA 污染
2. 用 DecontX 对模拟污染数据做去污染
3. 从表达量、污染比例估计、UMAP 和 Leiden 聚类几个角度评估校正效果

## 项目结构

建议把脚本放在项目根目录下的 scripts 文件夹中，因为代码里使用的是脚本所在位置的上一级目录作为项目根目录。

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

## 数据准备

原始 PBMC10k 数据放在下面这个目录：

```text
data/raw/pbmc10k
```

需要包含两个 10x h5 文件：

```text
10k_PBMC_3p_nextgem_Chromium_X_raw_feature_bc_matrix.h5
10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5
```

其中 filtered matrix 用来做基础 QC、聚类和模拟；raw matrix 主要用于检查文件是否能被 Scanpy 正常读取。

如果要运行 CellBender 对照部分，还需要提前准备：

```text
results/cellbender_final/pbmc10k_cellbender_filtered.h5
```

## 环境依赖

Python 部分主要用到：

```text
scanpy
numpy
pandas
scipy
scikit-learn
matplotlib
```

R 部分主要用到：

```text
Matrix
SingleCellExperiment
SummarizedExperiment
celda
```

DecontX 来自 celda 包。运行 R 脚本前需要确认这些 R 包已经安装好。

## 运行顺序

所有命令建议在项目根目录下运行。

### 1. 检查 PBMC10k 文件

```bash
python scripts/00_check_pbmc10k.py
```

这个脚本会检查 raw 和 filtered 两个 10x h5 文件是否存在，并确认 Scanpy 可以正常读取。

### 2. 生成干净基线数据

```bash
python scripts/01_qc_preprocess.py
```

这一步读取 filtered matrix，计算 QC 指标，过滤低质量细胞和基因，然后进行标准的 Scanpy 分析流程：

1. 过滤细胞和基因
2. 保存原始 counts 到 adata.layers["counts"]
3. normalize total 和 log1p
4. 选择 2000 个 highly variable genes
5. PCA、neighbors、UMAP
6. Leiden 聚类
7. marker gene 分析

主要输出：

```text
data/processed/pbmc10k_clean_baseline.h5ad
data/processed/pbmc10k_qc_metrics_before_filtering.csv
data/processed/pbmc10k_qc_metrics_after_filtering.csv
data/processed/pbmc10k_clean_baseline_marker_genes.csv
results/figures
```

后续模拟和 benchmark 都依赖 pbmc10k_clean_baseline.h5ad。

### 3. 模拟 ambient RNA 污染

```bash
python scripts/02_simulate_ambientRNA.py
```

这一步基于 clean baseline 里的原始 counts 来模拟 ambient RNA。脚本会挑选高表达、非线粒体、非核糖体基因作为 ambient profile。

模拟条件一共有六组：

```text
1gene_5pct
1gene_15pct
1gene_30pct
100genes_5pct
100genes_15pct
100genes_30pct
```

每个条件都会生成一个污染后的 h5ad 文件，并保存真实加入的 ambient counts，方便后面做 ground truth benchmark。

主要输出：

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

### 4. 先看不校正时的影响

```bash
python scripts/03_validate_simulation_no_correction.py
```

这一步把污染后的数据和 clean baseline 直接比较，用来确认模拟污染确实改变了表达和聚类结果。

比较内容包括：

1. observed contamination fraction
2. clean 和 contaminated 的 pseudobulk Pearson correlation
3. clean 和 contaminated 的 pseudobulk Spearman correlation
4. gene mean MAE 和 MSE
5. clean Leiden 和 contaminated Leiden 的 ARI、NMI
6. 每个污染条件下的 UMAP

主要输出：

```text
results/tables/pbmc10k_no_correction_simulation_benchmark.csv
results/figures/pbmc10k_observed_contamination_fraction.png
results/figures/pbmc10k_clean_vs_contaminated_pearson.png
results/figures/pbmc10k_clean_vs_contaminated_mse.png
results/figures/pbmc10k_clean_vs_contaminated_ari.png
results/figures/pbmc10k_clean_vs_contaminated_nmi.png
```

### 5. 导出 DecontX 输入文件

```bash
python scripts/04_export_decontx_inputs.py
```

DecontX 在 R 里运行，所以这里先把 Python 里的 h5ad 数据转换成 R 更方便读取的格式。

每个模拟条件会导出：

```text
counts_genes_by_cells.mtx
genes.csv
barcodes.csv
clusters.csv
```

clusters.csv 里的 cluster label 是用污染数据本身重新跑 Scanpy 聚类得到的，没有直接使用 clean baseline 的 ground truth cluster。

主要输出目录：

```text
data/decontx_input
```

### 6. 运行 DecontX

```bash
Rscript scripts/05_run_decontx_optimized.R
```

默认会对 data/decontx_input 下面的所有条件运行 DecontX。

也可以只跑指定条件，例如：

```bash
Rscript scripts/05_run_decontx_optimized.R 1gene_30pct 100genes_30pct
```

每个条件会输出 corrected counts 和每个 cell 的 DecontX contamination estimate。

主要输出目录：

```text
data/decontx_output
```

每个条件下主要包括：

```text
decontx_corrected_counts_genes_by_cells.mtx
decontx_cell_metadata.csv
genes.csv
barcodes.csv
```

### 7. 快速评估 DecontX 校正效果

```bash
python scripts/06_benchmark_decontx_quick.py
```

这一步比较 No correction 和 DecontX 两种结果相对于 clean baseline 的差异。

指标包括：

1. pseudobulk Pearson correlation
2. pseudobulk Spearman correlation
3. gene mean MAE
4. gene mean MSE
5. DecontX 估计的 mean 和 median contamination

主要输出：

```text
results/tables/pbmc10k_decontx_quick_benchmark.csv
results/figures/pbmc10k_decontx_quick_pearson.png
results/figures/pbmc10k_decontx_quick_mse.png
results/figures/pbmc10k_decontx_quick_mae.png
```

### 8. 画 DecontX 估计污染比例

```bash
python scripts/07_plot_decontx_contamination_estimate.py
```

这一步读取第 7 步生成的 benchmark 表，把设计的 contamination fraction 和 DecontX 估计值画在一起。

主要输出：

```text
results/figures/pbmc10k_decontx_estimated_vs_designed_contamination.png
```

### 9. 看 UMAP 和聚类是否恢复

```bash
python scripts/08_umap_clustering_impact.py
```

默认只分析两个高污染条件：

```text
1gene_30pct
100genes_30pct
```

也可以手动指定条件：

```bash
python scripts/08_umap_clustering_impact.py 1gene_5pct 1gene_15pct 1gene_30pct 100genes_5pct 100genes_15pct 100genes_30pct
```

这一步会对 contaminated data 和 DecontX corrected data 运行相同的 Scanpy downstream pipeline，然后和 clean baseline 的 Leiden 聚类做比较。

主要输出：

```text
results/tables/pbmc10k_decontx_umap_clustering_impact.csv
results/figures/pbmc10k_1gene_30pct_umap_by_clean_cluster.png
results/figures/pbmc10k_1gene_30pct_umap_by_own_leiden.png
results/figures/pbmc10k_100genes_30pct_umap_by_clean_cluster.png
results/figures/pbmc10k_100genes_30pct_umap_by_own_leiden.png
```

如果传入更多条件，会按条件生成对应的 UMAP 图。

### 10. CellBender 对照分析

```bash
python scripts/09_benchmark_cellbender.py
```

这一步是可选的，用来把 CellBender corrected matrix 和 Cell Ranger filtered baseline 做对照。

运行前需要有：

```text
results/cellbender_final/pbmc10k_cellbender_filtered.h5
```

分析内容包括：

1. Cell Ranger 和 CellBender 的 common cells、common genes
2. UMI counts 和 detected genes per cell 的变化
3. pseudobulk expression correlation
4. gene mean MAE 和 MSE
5. Leiden 聚类 ARI、NMI
6. UMAP 对比图
7. marker gene dotplot
8. 保存处理后的 AnnData 对象

主要输出：

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

## 一次性运行主流程

如果只跑 DecontX 主线，可以按下面顺序执行：

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

CellBender 对照单独运行：

```bash
python scripts/09_benchmark_cellbender.py
```

## 结果怎么看

主要表格在：

```text
results/tables
```

主要图片在：

```text
results/figures
```

最核心的结果可以先看这几个文件：

```text
results/tables/pbmc10k_no_correction_simulation_benchmark.csv
results/tables/pbmc10k_decontx_quick_benchmark.csv
results/tables/pbmc10k_decontx_umap_clustering_impact.csv
results/figures/pbmc10k_decontx_estimated_vs_designed_contamination.png
```

如果 DecontX 校正有效，通常会看到 DecontX 相比 No correction 更接近 clean baseline，例如 pseudobulk correlation 更高，gene-level error 更低，聚类结果和 clean baseline 的 ARI、NMI 更接近。

## 注意事项

1. Python 脚本默认从 scripts 文件夹往上找项目根目录，所以不要随意移动脚本位置。
2. R 脚本需要在项目根目录下运行，否则会找不到 data/decontx_input。
3. 01_qc_preprocess.py 必须先运行，因为后面的模拟需要 data/processed/pbmc10k_clean_baseline.h5ad。
4. 05_run_decontx_optimized.R 必须在 04_export_decontx_inputs.py 之后运行。
5. 07_plot_decontx_contamination_estimate.py 依赖 06_benchmark_decontx_quick.py 生成的 benchmark 表。
6. 09_benchmark_cellbender.py 不是 DecontX 主流程的一部分，只有在准备好 CellBender 输出后再运行。
