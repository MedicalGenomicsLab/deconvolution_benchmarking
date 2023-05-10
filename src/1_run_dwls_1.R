#### BUIDLING SIGNATURE MATRIX USING EITHER MAST OR SEURAT
# load packages
library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)

# Load dwls functions
source(
  "???/deconvolution_benchmarking_analysis/src/utils/1_DWLS_de_analysis.R"
)

temp_args <- commandArgs(trailingOnly = TRUE)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
cell_type <- temp_args[3]
output_dir <- temp_args[4]
sig_matrix_build_method <- temp_args[5]

# Load single-cell reference matrix
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
scRNA.ref.matrix <- as.matrix(scRNA.ref)

# Load single cell labels
sc.labels <- read.table(
  single_cell_labels_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)
sc_labels <- sc.labels$cell_labels

# Perform DE analysis using either Seurat or MAST for a specific cell type
if (sig_matrix_build_method == "seurat") {
  DEAnalysis(
    scdata = scRNA.ref.matrix,
    i = cell_type,
    id = sc_labels,
    path = output_dir
  )
} else if (sig_matrix_build_method == "mast") {
  DEAnalysisMAST(
    scdata = scRNA.ref.matrix,
    i = cell_type,
    id = sc_labels,
    path = output_dir
  )
} else {
  print("Unrecognised signature building method. Stop execution.")
  quit()
}
