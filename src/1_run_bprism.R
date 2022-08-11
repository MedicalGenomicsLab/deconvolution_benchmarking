# Load packages
library(TED)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
output_dir <- temp_args[4]
n_cpus <- as.numeric(temp_args[5])

# Load single-cell reference and bulk matrix
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
scRNA.ref.matrix <- as.matrix(scRNA.ref)

# Clean up ribosomal and gender genes
scRNA.ref.filtered <- cleanup.genes(
    ref.dat=scRNA.ref.matrix,
    species="hs", 
    gene.type=c("RB", "chrM", "chrX", "chrY"),
    input.type="scRNA", 
    exp.cells=5 
)

# Load single cell labels
sc.labels <- read.table(
  single_cell_labels_path,
  header = TRUE,
  sep="\t",
  row.names = 1,
  stringsAsFactors=FALSE
)

# Load bulk mixtures
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
# Tranpose and transform DataFrame to matrix
bulk.matrix <- as.matrix(t(bulk))

# Run BayesPrism
res <- run.Ted(
    ref.dat=scRNA.ref.filtered,
    X=bulk.matrix,
    cell.type.labels=sc.labels$cell_labels,
    tum.key="Cancer_Epithelial", 
    input.type="scRNA",
    n.cores=n_cpus,
    pdf.name=paste0(output_dir, "/results")
)

# Save output
write.csv(
  res$res$final.gibbs.theta, paste0(output_dir, "/results.csv")
)
save(
  res, file=paste0(output_dir, "/results.RData")
)
