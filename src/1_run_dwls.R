#load packages
library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)

# Load dwls functions
source(
  "/working/lab_nicw/khoaT/deep_tme/tme_profiling/tme_profiling_tools/dwls/R/functions.R"
)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
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

# Load single cell labels
sc.labels <- read.table(
  single_cell_labels_path,
  header = TRUE,
  sep="\t",
  row.names = 1,
  stringsAsFactors=FALSE
)
sc_labels <- sc.labels$cell_labels

# Load bulk mixtures
bulk <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
bulk.matrix <- as.matrix(bulk)

# Build signature from single-cell data
Signature<-buildSignatureMatrixUsingSeurat(
  scdata=scRNA.ref.matrix,
  id=sc_labels,
  path=output_dir,
  diff.cutoff=0.5,
  pval.cutoff=0.01
)

# DWLS only does deconvolution one bulk sample at a time. 
# Therefore we need to iterate over each column of the bulk matrix to estimate cell proportions
allCounts_DWLS<-NULL
for(j in 1:(dim(bulk.matrix)[2])){
  S<-Signature
  Bulk<-bulk.matrix[,j]
  names(Bulk)<-rownames(bulk.matrix)
  
  # Get gens intersection between bulk and signature matrix
  # Equivalent to trimData() in source code
  Genes<-intersect(rownames(S),names(Bulk))
  B<-Bulk[Genes]
  S<-S[Genes,]
  
  # Solve proportions using DWLS
  solDWLS<-solveDampenedWLS(S,B)
  
  # Append results
  allCounts_DWLS<-cbind(allCounts_DWLS, solDWLS)
}

# Retrieve bulk mixture ids and save predictions
colnames(allCounts_DWLS) <- colnames(bulk.matrix)
write.csv(
  allCounts_DWLS, paste0(output_dir, "/results.csv")
)
