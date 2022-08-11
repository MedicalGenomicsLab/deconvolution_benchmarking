library(scBio)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
single_cell_space_path <- temp_args[3]
bulk_mixture_path <- temp_args[4]
output_dir <- temp_args[5]
num_cores <- strtoi(temp_args[6])
neighborhood_size <- strtoi(temp_args[7])

# Load single-cell reference matrix
test.scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep="\t",
  row.names=1
)

# Load single cell labels
test.sc.labels <- read.table(
  single_cell_labels_path,
  header = TRUE,
  sep="\t",
  row.names = 1
)
sc_labels <- test.sc.labels$cell_labels

# Load single-cell space
test.sc.space <- read.table(
  single_cell_space_path,
  header=TRUE, 
  sep=",",
  row.names = 1
)
test.sc.space <- as.matrix(test.sc.space)

# Load bulk mixtures
test.bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)


# Run CPM in absolute mode (single-sample prediction)
res <- CPM(
  test.scRNA.ref, 
  sc_labels, 
  test.bulk, 
  test.sc.space,
  no_cores=num_cores,
  neighborhoodSize=neighborhood_size,
  quantifyTypes=TRUE,
  typeTransformation=TRUE
)

# Save cell-type predictions and also the whole RData object
write.csv(
  res$cellTypePredictions, paste0(output_dir, "/results.csv")
)
saveRDS(res, paste0(output_dir, "/results.rda"))
