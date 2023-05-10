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
  "???/src/utils/09_DWLS_de_analysis.R"
)
source(
  "???/src/utils/09_DWLS_deconvolution.R"
)

temp_args <- commandArgs(trailingOnly = TRUE)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
sig_matrix_build_method <- temp_args[4]
de_output_dir <- temp_args[5]
output_dir <- temp_args[6]

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

# Load bulk mixtures
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep = "\t",
  row.names = 1
)
bulk.matrix <- as.matrix(bulk)

# Build signature from single-cell data
if (sig_matrix_build_method == "seurat") {
  signature <- buildSignatureMatrixUsingSeurat(
    scdata = scRNA.ref.matrix,
    id = sc_labels,
    path = de_output_dir,
    diff.cutoff = 0.5,
    pval.cutoff = 0.01
  )
} else if (sig_matrix_build_method == "mast") {
  signature <- buildSignatureMatrixMAST(
    scdata = scRNA.ref.matrix,
    id = sc_labels,
    path = de_output_dir,
    diff.cutoff = 0.5,
    pval.cutoff = 0.01
  )
} else {
  print("Unrecognised signature building method. Stop execution.")
  quit()
}


# DWLS only does deconvolution one bulk sample at a time
# Therefore we need to iterate over each column of the bulk matrix to estimate cell proportions
allCounts_DWLS <- NULL
for (j in 1:(dim(bulk.matrix)[2])) {
  S <- signature
  Bulk <- bulk.matrix[, j]
  names(Bulk) <- rownames(bulk.matrix)

  # Get gens intersection between bulk and signature matrix
  # Equivalent to trimData() in source code
  Genes <- intersect(rownames(S), names(Bulk))
  B <- Bulk[Genes]
  S <- S[Genes, ]

  # Solve proportions using DWLS
  solDWLS <- solveDampenedWLS(S, B)

  # Append results
  allCounts_DWLS <- cbind(allCounts_DWLS, solDWLS)
}

# Retrieve bulk mixture ids and save predictions
colnames(allCounts_DWLS) <- colnames(bulk.matrix)
write.csv(
  allCounts_DWLS, paste0(output_dir, "/results.csv")
)
