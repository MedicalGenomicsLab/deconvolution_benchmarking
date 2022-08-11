library(jsonlite)
library(hspe)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
bulk_mixture_path <- temp_args[2]
pure_samples_path <- temp_args[3]
output_dir <- temp_args[4]

# Load single cell reference and bulk expressions
scRNA.ref <-read.table(
  single_cell_ref_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)

# Read pure samples
pure.samples <- fromJSON(pure_samples_path)

# Run hspe beautifully
# Here we don't specify method to select marker genes
# This means hspe will use its default method, i.e: top 10% of genes
hspe_out = hspe(
  Y=bulk, 
  reference=scRNA.ref,
  pure_samples=pure.samples,
  verbose=TRUE
)

# Save fraction estimations
write.csv(
  hspe_out$estimates, 
  paste0(output_dir, "/results.csv")
)
# Save entire results into .RData object
save(
  hspe_out, 
  file=paste0(output_dir, "/results.RData")
)

# Save cell markers
write(
  toJSON(hspe_out$markers), 
  file=paste0(output_dir, "/marker_genes.json")
)
