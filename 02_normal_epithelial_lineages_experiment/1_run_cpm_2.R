temp_args <- commandArgs(trailingOnly = T)
stage_1_output_path <- temp_args[1]
single_cell_ref_path <- temp_args[2]
bulk_mixture_path <- temp_args[3]
confidence_interval <- as.logical(temp_args[4])
output_dir <- temp_args[5]

source(
  "/working/lab_nicw/khoaT/deep_tme/tme_profiling/tme_benchmarking/deconvolution_benchmarking/09_CPM_functions.R"
)

# Load output from stage 1
load(stage_1_output_path)

# Confidence interval is already set to TRUE in stage 1
# Reset it to stage 2's input
calculateCI <- confidence_interval

# Load single cell reference and bulk expressions
# They need to be assigned to refference and Y specifically
refference <-read.table(
  single_cell_ref_path,
  header = TRUE,
  sep="\t",
  row.names=1
)
Y <- read.table(
  bulk_mixture_path,
  header = TRUE,
  sep="\t",
  row.names=1
)

##### Combining cell predictions #####
print("Combining CPM iterations")
predictedCells = matrix(0, nrow = dim(YReduced)[2], ncol = dim(refferenceSmaller)[2])
predictedCellsCounts = matrix(0, nrow = dim(YReduced)[2], ncol = dim(refferenceSmaller)[2])

for(resultMatrix in resultSmallMatrixes){
  completeResultMatrix = matrix(0, nrow = dim(resultMatrix)[2], ncol = dim(refferenceSmaller)[2])
  completeResultMatrix[,as.numeric(as.matrix(row.names(resultMatrix)))] = t(resultMatrix)
  predictedCells = predictedCells + completeResultMatrix
  predictedCellsCounts = predictedCellsCounts + abs(sign(completeResultMatrix))
}
predictedCellsFinal = predictedCells/predictedCellsCounts

##### Smoothing #####
print("Smoothing")
allClusterMeansMatrix = t(do.call(rbind,lapply(1:length(refferenceNames),function(cell){
  rowMeans(predictedCellsFinal[,chosenCellNeigList[[cell]]])
})))
colnames(allClusterMeansMatrix) = colnames(refference)
row.names(allClusterMeansMatrix) = colnames(Y)

cellTypeRes = NULL
seRes = NULL
confMatrix = NULL

#### Cell type prediction ####
if(quantifyTypes){
  print("Calculating cell type quantities")
  allClusterMeansMatrixForCellTypes = allClusterMeansMatrix
  if(typeTransformation){
    allClusterMeansMatrixForCellTypes = t(apply(t(allClusterMeansMatrixForCellTypes),2,function(x){
      x-min(x)
    }))
  }
  
  # cellTypeRes = do.call(cbind,lapply(unique(refferenceNames),function(currCluster){
  #   apply(allClusterMeansMatrixForCellTypes,1,function(x){
  #     #ks.test(x,x[currCluster==refferenceNames],alternative = "less")$p.value
  #     df = as.data.frame(cbind(sample(x,length(which(currCluster==refferenceNames))),x[currCluster==refferenceNames]))
  #     lm(V2~V1+0, data = df)$coefficients[1]
  #   })
  # }))
  
  cellTypeRes = do.call(cbind,lapply(unique(refferenceNames),function(currCluster){
    rowMeans(allClusterMeansMatrixForCellTypes[,currCluster==refferenceNames])
  }))
  colnames(cellTypeRes) = unique(refferenceNames)
  
  if(typeTransformation){
    cellTypeRes = t(apply(t(cellTypeRes),2,function(x){
      #x = (x-min(x))
      x/sum(x)
    })
    )
  }
}

#### Standard error prediction ####
if(calculateCI){
  print("Calculating the confidence interval matrix")
  
  resultOriginalSizeMatrixes = lapply(resultSmallMatrixes, function(resultSmallMatrix){
    completeResultMatrix = matrix(
      NA, 
      nrow = dim(resultSmallMatrix)[2], 
      ncol = dim(refferenceSmaller)[2]
    )
    completeResultMatrix[
      ,
      match(
        colnames(allClusterMeansMatrix)[
          as.numeric(as.matrix(row.names(resultSmallMatrix)))
        ],
        colnames(refferenceSmaller)
      )
    ] = t(resultSmallMatrix)
    completeResultMatrix
  })
  
  seRes <- do.call(rbind,lapply(colnames(YReduced), function(sample){
    sampleMatrix = do.call(rbind, lapply(resultOriginalSizeMatrixes,function(currRes){
      currRes[which(colnames(YReduced)==sample),]
    }))
    apply(sampleMatrix,2,function(x){
      sd(x[!is.na(x)])/sqrt(length(which(!is.na(x))))
    })
  }))
  
  seResNorm = t(do.call(rbind,lapply(1:length(refferenceNames),function(cell){
    rowMeans(seRes[,chosenCellNeigList[[cell]]])
  })))
  
  confMatrix = matrix(paste(allClusterMeansMatrix-1.96*seResNorm,allClusterMeansMatrix+1.96*seResNorm,sep = " <-> "),ncol = dim(allClusterMeansMatrix)[2])
  
  colnames(seRes) = colnames(confMatrix) = colnames(refference)
  row.names(seRes) = row.names(confMatrix) = colnames(Y)
}

print("Done")
results <-list(predictions = allClusterMeansMatrix, cellTypePredictions = cellTypeRes, sePredictions = seRes, confMatrix = confMatrix)

# Save cell-type predictions into csv and results into an RData object
write.csv(
  results$cellTypePredictions, 
  paste0(output_dir, "/results.csv")
)
saveRDS(
  results, 
  file=paste0(output_dir, "/stage_2_output.RData")
)
