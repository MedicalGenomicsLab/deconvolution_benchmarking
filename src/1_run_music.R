library(MuSiC)
library(Biobase)
library(xbioc)
library(jsonlite)
library(testit)

temp_args <- commandArgs(trailingOnly = TRUE)
single_cell_ref_path <- temp_args[1]
single_cell_phenotypes_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
use_marker_genes <- as.logical(temp_args[4])
marker_genes_path <- temp_args[5]
use_cell_subtypes <- as.logical(temp_args[6])
cell_subtypes_path <- temp_args[7]
output_dir <- temp_args[8]

# Load single-cell reference matrix
print("Loading single-cell reference")
scRNA.ref <- read.table(
  single_cell_ref_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)
scRNA.ref.matrix <- as.matrix(scRNA.ref)

# Load single cell phenotypes
sc.pheno <- read.table(
  single_cell_phenotypes_path,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  as.is = c("cellTypeID", "sampleID")
)
sc.pheno.meta <- data.frame(
  labelDescription = c(
    "SubjectName", "cellType", "cellTypeID", "sampleID"
  ),
  row.names = c(
    "Subject Name", "cell Type", "cell Type ID", "sample ID"
  )
)

sc.pdata <- new(
  "AnnotatedDataFrame",
  data = sc.pheno, varMetadata = sc.pheno.meta
)
sc.eset <- Biobase::ExpressionSet(
  assayData = scRNA.ref.matrix, phenoData = sc.pdata
)

# Load bulk mixtures and put it in an ExpressionSet object
print("Loading bulk mixtures")
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep = "\t",
  row.names = 1
)
bulk.matrix <- as.matrix(bulk)
bulk.eset <- Biobase::ExpressionSet(assayData = bulk.matrix)

# If we're using cell subtypes, we'll run MuSiC with music_prop.cluster()
if (use_cell_subtypes == TRUE) {
  print("Running with cell-subtypes")

  # use_marker_genes has to be TRUE. If not, terminate execution
  # marker_genes_path and cell_subtypes_path also have to be specified
  assert("Marker genes are needed for music_prop.cluster()", {
    (use_marker_genes <- TRUE)
    (marker_genes_path != "None")
    (cell_subtypes_path != "None")
  })
  # Load marker genes
  marker.genes <- fromJSON(marker_genes_path)

  # Load cluster types
  clusters.type <- fromJSON(cell_subtypes_path)

  # Create ClusterType attribute in single-cell ExpressionSet
  cl.type <- as.character(sc.eset$cellType)
  for (cl in 1:length(clusters.type)) {
    cl.type[cl.type %in% clusters.type[[cl]]] <- names(clusters.type)[cl]
  }
  pData(sc.eset)$clusterType <- factor(
    cl.type,
    levels = c(names(clusters.type))
  )

  # Run MuSiC with pre-defined clustering
  music.res <- music_prop.cluster(
    bulk.eset = bulk.eset,
    sc.eset = sc.eset,
    group.markers = marker.genes,
    groups = "clusterType",
    clusters = "cellType",
    clusters.type = clusters.type,
    samples = "sampleID",
  )

  # Save cell-type predictions and also the whole RData object
  write.csv(
    music.res$Est.prop.weighted.cluster, paste0(output_dir, "/results.csv")
  )
  saveRDS(music.res, paste0(output_dir, "/results.rds"))
} else {
  print("Running without cell-subtypes")

  # Check if marker genes was specified
  if (use_marker_genes == TRUE) {
    print("Running with marker genes")

    # marker_genes_path has to be specified
    assert("Marker genes are needed for music_prop.cluster()", {
      (marker_genes_path != "None")
    })

    # Load marker genes
    marker.genes <- fromJSON(marker_genes_path)
  } else {
    print("Running without marker genes")

    # Assign NULL to marker genes
    marker.genes <- NULL
  }

  # Run MuSiC with normal mode
  music.res <- music_prop(
    bulk.eset = bulk.eset,
    sc.eset = sc.eset,
    clusters = "cellType",
    samples = "sampleID",
    select.ct = as.vector(unique(sc.pheno$cellType)),
    markers = marker.genes
  )

  # Save cell-type predictions and also the whole RData object
  write.csv(
    music.res$Est.prop.weighted, paste0(output_dir, "/results.csv")
  )
  saveRDS(music.res, paste0(output_dir, "/results.rds"))
}
