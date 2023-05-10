# Load packages
library(BayesPrism)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
output_dir <- temp_args[4]
n_cpus <- as.numeric(temp_args[5])
use_protein_coding_genes <- as.logical(temp_args[6])
use_marker_genes <- as.logical(temp_args[7])
use_cell_state_labels <- as.logical(temp_args[8])
cancer_label <- temp_args[9]

# Load single-cell reference and bulk matrix
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
scRNA.ref.matrix <- as.matrix(scRNA.ref)

# Clean up ribosomal and gender genes
scRNA.ref.filtered <- cleanup.genes(
  input = scRNA.ref.matrix,
  input.type = "count.matrix",
  species = "hs",
  gene.group = c("Rb", "Mrp", "other_Rb", "chrM", "MALAT1", "chrX", "chrY"),
  exp.cells = 5
)

# Load single cell labels
sc.labels <- read.table(
  single_cell_labels_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)
# Get cell-type and cell-state labels
cell.type.labels <- sc.labels$cell_type_labels
if (use_cell_state_labels == TRUE) {
  cell.state.labels <- sc.labels$cell_state_labels
} else {
  cell.state.labels <- NULL
}

# Load bulk mixtures
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep = "\t",
  row.names = 1
)
# Tranpose and transform DataFrame to matrix
bulk.matrix <- as.matrix(t(bulk))


# If use_protein_coding_genes is TRUE
# Filter our non-protein coding genes
if (use_protein_coding_genes == TRUE) {
  scRNA.ref.filtered <- select.gene.type(
    scRNA.ref.filtered,
    gene.type = "protein_coding"
  )
}
# If use_marker_genes is TRUE
# Conduct t-test and filter marker genes based on results
if (use_marker_genes == TRUE) {
  diff.exp.stat <- get.exp.stat(
    sc.dat = scRNA.ref.matrix[, colSums(scRNA.ref.matrix > 0) > 3], # filter genes to reduce memory use
    cell.type.labels = cell.type.labels,
    cell.state.labels = cell.state.labels,
    psuedo.count = 0.1, # a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
    cell.count.cutoff = 50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
    n.cores = n_cpus # number of threads
  )
  scRNA.ref.filtered <- select.marker(
    sc.dat = scRNA.ref.filtered,
    stat = diff.exp.stat,
    pval.max = 0.01,
    lfc.min = 0.1
  )
}


# Run BayesPrism
myPrism <- new.prism(
  reference = scRNA.ref.filtered,
  mixture = bulk.matrix,
  input.type = "count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key = cancer_label,
  outlier.cut = 0.01,
  outlier.fraction = 0.1,
)
bp.res <- run.prism(prism = myPrism, n.cores = n_cpus)


# Save cell-type fractions and also whole BayesPrism output object
theta <- get.fraction(
  bp = bp.res,
  which.theta = "final",
  state.or.type = "type"
)
write.csv(
  theta, paste0(output_dir, "/results.csv")
)
save(
  bp.res,
  file = paste0(output_dir, "/results.RData")
)
