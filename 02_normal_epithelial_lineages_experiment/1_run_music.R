library(MuSiC)
library(Biobase)
library(xbioc)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_phenotypes_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
output_dir <- temp_args[4]


# Load single-cell reference matrix
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
scRNA.ref.matrix <- as.matrix(scRNA.ref)

# Load single cell phenotypes
sc.pheno <- read.table(
  single_cell_phenotypes_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
sc.pheno.meta <- data.frame(
  labelDescription=c(
    "SubjectName", "cellType", "cellTypeID", "sampleID"
  ),
  row.names=c(
    "Subject Name", "cell Type", "cell Type ID", "sample ID"
  )
)

sc.pdata <- new(
  "AnnotatedDataFrame", data=sc.pheno, varMetadata=sc.pheno.meta
)
sc.eset <- Biobase::ExpressionSet(
  assayData=scRNA.ref.matrix, phenoData=sc.pdata
)

# Load bulk mixtures and put it in an ExpressionSet object
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
bulk.matrix <- as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

# Run MuSiC beautifully
music.res <- music_prop(
  bulk.eset = bulk.eset, 
  sc.eset = sc.eset,
  clusters = 'cellType', 
  samples = 'sampleID', 
  select.ct = as.vector(unique(sc.pheno$cellType))
)

# Save cell-type predictions and also the whole RData object
write.csv(
  music.res$Est.prop.weighted, paste0(output_dir, "/results.csv")
)
saveRDS(music.res, paste0(output_dir, "/results.rda"))
