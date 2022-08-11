library(EPIC)

temp_args <- commandArgs(trailingOnly = T)
ref_path <- temp_args[1]
bulk_mixture_path <- temp_args[2]
sig_matrix <- temp_args[3]
sig_genes_path <- temp_args[4]
output_dir <- temp_args[5]

if (sig_matrix=="default") {
  # Get reference profiles
  genes.ref <- load(ref_path)

  test.bulk <- read.table(
    bulk_mixture_path,
    header = TRUE,
    sep="\t",
    row.names=1
  )

  results <- EPIC(bulk = test.bulk, reference = genes.ref)

} else if (sig_matrix=="cbx") {
  # Get reference profiles
  ref.profiles <- read.table(
    ref_path,
    header = TRUE,
    sep="\t",
    row.names=1
  )
  ref.profiles.matrix <- as.matrix(ref.profiles)

  # Get signature genes
  sig.genes.matrix <- read.table(
    sig_genes_path,
    header = TRUE,
    sep="\t",
    row.names=1
  )
  sig.genes <- sig.genes.matrix$gene_symbol

  # Put reference profiles and signature genes in a list
  ref.profiles.list <- list(
    refProfiles=ref.profiles.matrix, sigGenes=sig.genes
  )

# Load bulk mixtures
  test.bulk <- read.table(
    bulk_mixture_path,
    header = TRUE,
    sep="\t",
    row.names=1
  )
  test.bulk.matrix <- as.matrix(test.bulk)

  results <- EPIC(bulk = test.bulk.matrix, reference = ref.profiles.list)

} else {
  print("Unrecognized signature matrix")
}

# Save cell-type predictions and also the whole RData object
write.csv(
  results$cellFractions, paste0(output_dir, "/results.csv")
)
saveRDS(results, paste0(output_dir, "/results.rda"))
