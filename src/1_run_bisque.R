library(BisqueRNA)
library(Biobase)
library(testit)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_phenotypes_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
use_marker_genes <- as.logical(temp_args[4])
marker_genes_path <- temp_args[5]
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

# Load single cell phenotypes
sc.pheno.matrix <- read.table(
  single_cell_phenotypes_path,
  header = TRUE,
  sep = "\t"
)
sample.ids <- sc.pheno.matrix$cell_ids
individual.labels <- sc.pheno.matrix$patient_ids
cell.type.labels <- sc.pheno.matrix$cell_labels

# Combine phenotype and single cell counts into an ExpressionSet object
sc.pheno <- data.frame(
  check.names = FALSE,
  check.rows = FALSE,
  stringsAsFactors = FALSE,
  row.names = sample.ids,
  SubjectName = individual.labels,
  cellType = cell.type.labels
)
sc.meta <- data.frame(
  labelDescription = c("SubjectName", "cellType"),
  row.names = c("SubjectName", "cellType")
)
sc.pdata <- new(
  "AnnotatedDataFrame",
  data = sc.pheno, varMetadata = sc.meta
)
sc.eset <- Biobase::ExpressionSet(
  assayData = scRNA.ref.matrix, phenoData = sc.pdata
)

# Load bulk mixtures and put it in an ExpressionSet object
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep = "\t",
  row.names = 1
)
bulk.matrix <- as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

# Check if marker genes was specified
if (use_marker_genes == TRUE) {
  print("Running with marker genes")

  # marker_genes_path has to be specified
  assert("Marker genes are needed for music_prop.cluster()", {
    (marker_genes_path != "None")
  })

  marker.genes <- read.table(
    marker_genes_path,
    header = TRUE,
    sep = "\t"
  )
} else {
  marker.genes <- NULL
}

# Run bisque beautifully
res <- BisqueRNA::ReferenceBasedDecomposition(
  bulk.eset, sc.eset,
  markers = marker.genes, use.overlap = FALSE
)

# Save cell-type predictions and also the whole RData object
write.csv(
  res$bulk.props, paste0(output_dir, "/results.csv")
)
saveRDS(res, paste0(output_dir, "/results.rda"))
