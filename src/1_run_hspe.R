library(jsonlite)
library(hspe)
library(testit)

temp_args <- commandArgs(trailingOnly = TRUE)
single_cell_ref_path <- temp_args[1]
bulk_mixture_path <- temp_args[2]
pure_samples_path <- temp_args[3]
use_marker_genes <- as.logical(temp_args[4])
marker_genes_path <- temp_args[5]
output_dir <- temp_args[6]

# Load single cell reference and bulk expressions
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Read pure samples
pure.samples <- fromJSON(pure_samples_path)

# Check if marker genes was specified
if (use_marker_genes == TRUE) {
  print("Running with marker genes")

  # marker_genes_path has to be specified
  assert("Marker genes are needed when use_marker_genes is TRUE", {
    (marker_genes_path != "None")
  })

  # Load marker genes
  marker.genes <- fromJSON(marker_genes_path)
} else {
  print("Running without marker genes")

  # Assign NULL to marker genes
  marker.genes <- NULL
}

# Run hspe beautifully
# Here we don't specify method to select marker genes
# This means hspe will use its default method, i.e: top 10% of genes
hspe_out <- hspe(
  Y = bulk,
  reference = scRNA.ref,
  pure_samples = pure.samples,
  markers = marker.genes,
  verbose = TRUE
)

# Save fraction estimations
write.csv(
  hspe_out$estimates,
  paste0(output_dir, "/results.csv")
)
# Save entire results into .RData object
save(
  hspe_out,
  file = paste0(output_dir, "/results.RData")
)

# Save cell markers
write(
  toJSON(hspe_out$markers),
  file = paste0(output_dir, "/marker_genes.json")
)
