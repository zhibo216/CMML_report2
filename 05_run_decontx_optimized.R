suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(celda)
})

# Important:
# Run this script from the project root:
# D:/CMML2/project8_ambientRNA

args <- commandArgs(trailingOnly = TRUE)

base_dir <- normalizePath(getwd(), winslash = "/")
input_base <- file.path(base_dir, "data", "decontx_input")
output_base <- file.path(base_dir, "data", "decontx_output")

if (!dir.exists(input_base)) {
  stop("Cannot find input folder: ", input_base,
       "\nPlease run this from the project root folder.")
}

dir.create(output_base, recursive = TRUE, showWarnings = FALSE)

if (length(args) >= 1) {
  conditions <- args
} else {
  conditions <- list.dirs(input_base, recursive = FALSE, full.names = FALSE)
}

cat("Conditions to run:\n")
print(conditions)

for (condition in conditions) {
  cat("\n==============================\n")
  cat("Running DecontX for:", condition, "\n")
  cat("==============================\n")

  in_dir <- file.path(input_base, condition)
  out_dir <- file.path(output_base, condition)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  counts_path <- file.path(in_dir, "counts_genes_by_cells.mtx")
  genes_path <- file.path(in_dir, "genes.csv")
  barcodes_path <- file.path(in_dir, "barcodes.csv")
  clusters_path <- file.path(in_dir, "clusters.csv")

  if (!file.exists(counts_path)) {
    stop("Missing counts file: ", counts_path)
  }
  if (!file.exists(genes_path)) {
    stop("Missing genes file: ", genes_path)
  }
  if (!file.exists(barcodes_path)) {
    stop("Missing barcodes file: ", barcodes_path)
  }
  if (!file.exists(clusters_path)) {
    stop("Missing clusters file: ", clusters_path)
  }

  counts <- readMM(counts_path)
  genes <- read.csv(genes_path, stringsAsFactors = FALSE)
  barcodes <- read.csv(barcodes_path, stringsAsFactors = FALSE)
  clusters <- read.csv(clusters_path, stringsAsFactors = FALSE)

  rownames(counts) <- genes$gene
  colnames(counts) <- barcodes$barcode

  # Match cluster order to count matrix columns
  clusters <- clusters[match(colnames(counts), clusters$barcode), ]

  if (any(is.na(clusters$decontx_input_cluster))) {
    stop("Some cells do not have DecontX input cluster labels.")
  }

  z <- as.character(clusters$decontx_input_cluster)

  cat("Input matrix dimensions:\n")
  cat("Genes:", nrow(counts), "\n")
  cat("Cells:", ncol(counts), "\n")
  cat("Number of input clusters:", length(unique(z)), "\n")

  sce <- SingleCellExperiment(
    assays = list(counts = counts)
  )

  set.seed(12345)

  # DecontX input:
  # counts = observed contaminated counts
  # z      = preliminary cell population labels from contaminated data
  sce <- decontX(
    sce,
    z = z,
    seed = 12345,
    verbose = TRUE
  )

  corrected <- decontXcounts(sce)

  if (!inherits(corrected, "sparseMatrix")) {
    corrected <- Matrix(corrected, sparse = TRUE)
  }

  corrected_path <- file.path(
    out_dir,
    "decontx_corrected_counts_genes_by_cells.mtx"
  )

  writeMM(corrected, corrected_path)

  coldata <- as.data.frame(colData(sce))

  contamination_col <- grep(
    "contamination",
    colnames(coldata),
    value = TRUE
  )[1]

  if (is.na(contamination_col)) {
    contamination_values <- rep(NA, ncol(sce))
  } else {
    contamination_values <- coldata[[contamination_col]]
  }

  metadata <- data.frame(
    barcode = colnames(sce),
    decontx_input_cluster = z,
    decontx_contamination = contamination_values
  )

  write.csv(
    metadata,
    file.path(out_dir, "decontx_cell_metadata.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(gene = rownames(sce)),
    file.path(out_dir, "genes.csv"),
    row.names = FALSE
  )

  write.csv(
    data.frame(barcode = colnames(sce)),
    file.path(out_dir, "barcodes.csv"),
    row.names = FALSE
  )

  cat("Saved corrected counts to:\n")
  cat(corrected_path, "\n")

  cat("Saved metadata to:\n")
  cat(file.path(out_dir, "decontx_cell_metadata.csv"), "\n")
}

cat("\nAll DecontX jobs finished.\n")
