library(scBio)

temp_args <- commandArgs(trailingOnly = T)
single_cell_ref_path <- temp_args[1]
single_cell_labels_path <- temp_args[2]
single_cell_space_path <- temp_args[3]
bulk_mixture_path <- temp_args[4]
output_dir <- temp_args[5]
num_cores <- strtoi(temp_args[6])
neighborhood_size <- strtoi(temp_args[7])

source(
  "/working/lab_nicw/khoaT/deep_tme/tme_profiling/tme_benchmarking/deconvolution_benchmarking/src/1_CPM_functions.R"
)

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


# Run first stage of CPM
SCData <- test.scRNA.ref
SCLabels <- sc_labels
BulkData <- test.bulk
cellSpace <- test.sc.space
no_cores <- num_cores
neighborhoodSize <- neighborhood_size
modelSize <- 10
minSelection <- 5
quantifyTypes <- TRUE
typeTransformation <- TRUE
calculateCI <- TRUE


genePercents = 0.4
if(min(table(SCLabels))<neighborhoodSize){
  neighborhoodSize = min(table(SCLabels))
  print(paste("Neighborhood size was switched to:",neighborhoodSize,sep=" "))
}
if(quantifyTypes){
  modelSize = length(unique(SCLabels))
  print(paste("Model size was switched to:",modelSize,sep=" "))
}else if(length(SCLabels)<modelSize){
  modelSize = length(SCLabels)
  print(paste("Model size was switched to:",modelSize,sep=" "))
}
if(!is.null(SCData) & !is.null(SCLabels) & !is.null(BulkData) & !is.null(cellSpace)){
  print("Selecting cells for each iteration")
}
cellSelection = choseCellsForRuns(cellSpace, SCLabels, modelSize, minSelection,neighborhoodSize)
numOfRunsToUse = cellSelection$numOfRuns
print(paste("Number of iteration:",numOfRunsToUse,sep=" "))
cellSelectionList = cellSelection$chosenCellList
cellNeigSelectionList = cellSelection$chosenNeigList

print("Running CPM, this may take a few minutes")


# Re-assign variables to match content of the CPMMAIN() function
refference <- SCData
refferenceNames <- SCLabels
Y <- BulkData
chosenCellList <- cellSelectionList
chosenCellNeigList <- cellNeigSelectionList
numOfRuns <- numOfRunsToUse

### RUN CPMMAIN() FUNCTION ###
YReduced = Y[row.names(Y) %in% row.names(refference), , drop = FALSE]

##### Revome genes low in reference data  #####
geneVarianceRef = apply(refference,1,function(gene){checkVariableGenes(as.numeric(as.matrix(gene)),0.1)})
geneVarianceFinalRef = sort(geneVarianceRef[geneVarianceRef>0],decreasing = T)
mutualGenes = names(geneVarianceFinalRef)[names(geneVarianceFinalRef) %in% row.names(YReduced)]

YReduced = YReduced[mutualGenes, , drop = FALSE]
refferenceSmaller = refference[mutualGenes,]

##### Main algorithm runs #####
if(is.null(no_cores)){
  no_cores = max(1, parallel::detectCores() - 1)
}
cl<-parallel::makeCluster(no_cores)
parallel::clusterExport(
  cl=cl, 
  varlist=c(
    "refferenceNames", "refferenceSmaller", "YReduced","neighborhoodSize","modelSize",
    "createSpecificRef","GeneBasedAnova", "chosenCellList", "chosenCellNeigList" ,"createNoCellDupeReference",
    "selectGenesUsingKappa", "runLibLinear", "genePercents"
  ),
  envir=environment()
)
doSNOW::registerDoSNOW(cl)
pb <- utils::txtProgressBar(min = 1, max = numOfRuns, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
#t = proc.time()
`%dopar2%` <- foreach::`%dopar%`
runNumber = NULL
resultSmallMatrixes <- foreach::foreach(runNumber = 1:numOfRuns, .options.snow = opts) %dopar2% {
  #for(runNumber in 1:numOfRuns) {
  print(runNumber)
  completeSpecificRefBefore = createSpecificRef(refferenceSmaller, modelSize, neighborhoodSize, genePercents, chosenCellList[[runNumber]], chosenCellNeigList)
  completeSpecificRef = completeSpecificRefBefore$ref
  
  #list of genes between clusters
  clusterNamesVector = rep("", dim(completeSpecificRef)[2])
  for(cluster in unique(refferenceNames)){
    selectedCellsForCluster = completeSpecificRefBefore$chosenCells[which(refferenceNames[completeSpecificRefBefore$chosenCells] == cluster)]
    selectedNamesForCluster = paste(colnames(refferenceSmaller)[selectedCellsForCluster],selectedCellsForCluster,sep="_")
    clusterNamesVector[!is.na(match(colnames(completeSpecificRef),selectedNamesForCluster))] = cluster
  }
  
  allGenes = c()
  
  #list of genes inside clusters
  
  if(quantifyTypes){
    allGenesList = GeneBasedAnova(completeSpecificRef)
  }else{
    allGenesList = lapply(unique(refferenceNames), function(cluster){
      specificClusterRef = completeSpecificRef[,which(clusterNamesVector == cluster)]
      colnames(specificClusterRef) = colnames(completeSpecificRef)[clusterNamesVector == cluster]
      GeneBasedAnova(specificClusterRef)
    })
  }
  
  for (list in allGenesList){
    allGenes = c(allGenes, list)
  }
  
  specificRefNoDupes = createNoCellDupeReference(completeSpecificRef)
  
  results = selectGenesUsingKappa(specificRefNoDupes, allGenes)
  X = results$reference
  X = refferenceSmaller[row.names(X),completeSpecificRefBefore$chosenCells]
  X = X[rowSums(X)!=0,]
  
  YRefinedReduced = YReduced[row.names(YReduced) %in% row.names(X),]
  PBSReductionData = YRefinedReduced
  
  setTxtProgressBar(pb, runNumber)
  
  resMatrix = t(runLibLinear(X, PBSReductionData))
  row.names(resMatrix) = chosenCellList[[runNumber]]
  resMatrix
}
#print(proc.time() - t)
parallel::stopCluster(cl)
close(pb)

save(
  chosenCellNeigList,
  Y,
  YReduced, 
  refference,
  refferenceNames,
  refferenceSmaller, 
  resultSmallMatrixes, 
  quantifyTypes,
  typeTransformation,
  calculateCI,
  file=paste0(output_dir, "/stage_1_output.rda")
)
