#dependencies
library(foreach)
library(doSNOW)
library(parallel)
library(raster)
library(fields)
library(limma)
library(LiblineaR)
library(sp)

########## create sub-reference for the deconvolution
#' @keywords internal
createSpecificRef <- function(currRefference, modelSize, neighborhoodSize, genePercents, chosenCells, chosenNeigCells){
  specificRef = do.call(cbind,lapply(chosenCells,function(chosenCell){
    currRefference[,chosenNeigCells[[chosenCell]]]
  }))
  colnames(specificRef) = unlist(lapply(chosenCells,function(cell){
    rep(paste(colnames(currRefference)[cell],cell,sep="_"),neighborhoodSize)
  }))
  row.names(specificRef) = row.names(currRefference)
  specificRef = specificRef[sample(1:dim(currRefference)[1],round(genePercents*dim(currRefference)[1])),]
  list(ref = specificRef, chosenCells = chosenCells)
}

########## Remove cell duplications from reference
#' @keywords internal
createNoCellDupeReference <- function(refference){
  currCellNames = colnames(refference)
  refferenceNoDups = as.data.frame(matrix(0,nrow = dim(refference)[1], ncol=length(unique(currCellNames))))
  for (i in 1:length(unique(currCellNames))){
    if (length(which(currCellNames == unique(currCellNames)[i]))!=1){
      refferenceNoDups[,i] = rowMeans(refference[,which(currCellNames == unique(currCellNames)[i])])
    }else{
      refferenceNoDups[,i] = refference[,which(currCellNames == unique(currCellNames)[i])]
    }
  }
  row.names(refferenceNoDups) = row.names(refference)
  colnames(refferenceNoDups) = unique(currCellNames)
  return(refferenceNoDups)
}

########## Score genes using ANOVA
#' @keywords internal
GeneBasedAnova <- function(specificRefference, nonZeroRatio = NULL){
  cellNamesForAnova = colnames(specificRefference)
  #genes_to_take = row.names(specificRefference)[!(grepl("Rpl", row.names(specificRefference)) | grepl("Rps", row.names(specificRefference)) | grepl("mt-", row.names(specificRefference)) | grepl("Gm[0-9]", row.names(specificRefference)))]
  genes_to_take = row.names(specificRefference)
  genes_to_take = names(which(rowMeans(specificRefference[genes_to_take,])!=0))
  if(!is.null(nonZeroRatio)){
    genes_to_take = unlist(apply(as.data.frame(genes_to_take),1,function(gene){
      if((length(which(as.numeric(specificRefference[gene,which(cellNamesForAnova==1)])!=0))/length(as.numeric(specificRefference[gene,which(cellNamesForAnova==1)])))>nonZeroRatio){
        gene
      }else{
        NULL
      }
    }))
  }
  
  dat = cbind(rep(0,length(genes_to_take)),specificRefference[genes_to_take,])
  group = c("",cellNamesForAnova)
  #dat = specificRefference[genes_to_take,]
  #group = c(cellNamesForAnova)
  #print(cellNamesForAnova)
  dmat <- stats::model.matrix(~ group)
  # fit the ANOVA model
  fit <- limma::lmFit(dat, dmat)
  fit = fit[,-1]
  fit <- limma::eBayes(fit)
  fitF = fit$F
  res = as.data.frame(cbind(gsub("group","",colnames(fit$coefficients)[apply(fit$coefficients,1,function(x){order(x,decreasing = T)[1]})]),fitF))
  colnames(res) = c("group", "score")
  
  listOfGenes = apply(as.data.frame(unique(cellNamesForAnova)),1,function(cellGroup){
    selectedIndexes = which(as.character(res$group)==as.character(cellGroup))
    (genes_to_take[selectedIndexes])[order(res$score[selectedIndexes],decreasing = T)]
  })
  return(listOfGenes)
}

########## Select final gene sets for reference using the kappa function
#' @keywords internal
selectGenesUsingKappa <- function(refferenceNoDups, allGenes){
  bestKappa = Inf
  bestG = 0
  mul = 1
  maxNumberOfGenesPerCell = 50
  bestGenes = c()
  indexRange = 2:maxNumberOfGenesPerCell
  for (i in indexRange){
    selectedGenes = unique(as.character(unlist(lapply(1:length(allGenes), function(listIndex){
      unlist(allGenes[listIndex])[which(!is.na(unlist(allGenes[listIndex])[1:as.numeric(i*mul)]))]
    }))))
    currRefferenceNoDups = refferenceNoDups[match(selectedGenes, row.names(refferenceNoDups)),]
    newKappa = kappa(currRefferenceNoDups)
    if (newKappa<bestKappa){
      bestKappa = newKappa
      bestG = i
      bestGenes = unique(selectedGenes)
    }
  }
  finalRefference = refferenceNoDups[which(row.names(refferenceNoDups) %in% bestGenes),]
  
  return(list(reference = finalRefference, G = bestG, kappa = bestKappa))
}

########## Check the variance of none zero genes
#' @keywords internal
checkVariableGenes = function(a, ratio) {
  count_nonZeros = length(which(a > min(a)))
  if (count_nonZeros/length(a) > ratio) {
    var(a)/ mean(a)
  } else {
    0
  }
}

########## Run svr based deconvolution
#' @keywords internal
runLibLinear = function(ref_matrix, sample_matrix){
  X <- data.matrix(ref_matrix)
  X <- X[apply(X,1,sd)!=0,]
  Y <- data.matrix(sample_matrix)
  Y = Y[match(row.names(X),row.names(Y)),]
  
  #intersect genes
  Y <- Y[row.names(Y) %in% row.names(X),]
  X <- X[row.names(X) %in% row.names(Y),]
  
  #standardize sig matrix
  #X <- (X - mean(X)) / sd(as.vector(X))
  X <- t(apply(X,1,function(rowX){
    (rowX - mean(rowX)) / sd(as.vector(rowX))
  }))
  
  C = LiblineaR::heuristicC(X)
  
  #iterate through mixtures
  predictionMatrix = do.call(rbind,lapply(1:dim(Y)[2],function(index){
    y <- Y[,index]
    
    #standardize mixture
    y <- (y - mean(y)) / sd(y)
    
    #run SVR core algorithm
    (LiblineaR::LiblineaR(data = X, target = y, type = 11, cost = C)$W)[1:dim(X)[2]]
  }))
  colnames(predictionMatrix) = colnames(X)
  row.names(predictionMatrix) = colnames(Y)
  predictionMatrix
}

########## Select cells for each run and caculate the desired number of runs
#' @keywords internal
choseCellsForRuns = function(XY, refNames, modelSize, minSelection, neighborhoodSize){
  #### cell selection
  k = floor(modelSize/length(unique(refNames)))
  if(k==0){
    k=1
  }
  
  minValueToReduceTo = 10^-10
  
  initialGrids = lapply(unique(refNames), function(currCluster){
    clusterIndexes = which(refNames==currCluster)
    nbins = max(k,length(clusterIndexes)/neighborhoodSize)
    if(is.null(dim(XY))){
      currXY = XY[clusterIndexes]
      breaks = seq(min(currXY)-10^-7,max(currXY)+10^-7, (max(currXY)-min(currXY)+2*10^-7)/ceiling(nbins))
      grid <- rep(NA,ceiling(nbins))
      cellLocationOnGrid = rep(NA,length(currXY))
      for(currBreakIndex in 2:length(breaks)){
        cellLocationOnGrid[which(currXY>breaks[currBreakIndex-1] & currXY<breaks[currBreakIndex])] = currBreakIndex-1
      }
      tab <- table(cellLocationOnGrid)
      grid[as.numeric(names(tab))] <- tab
    }else{
      currXY = XY[clusterIndexes,]
      ch <- grDevices::chull(currXY)
      coords <- currXY[c(ch, ch[1]), ]
      # ch = geometry::convhulln(currXY[,1:3])
      # coords = t(apply(ch,1,function(currCh){
      #   unlist(lapply(1:length(currCh),function(currIndex){
      #     currXY[currCh[currIndex],currIndex]
      #   }))
      # }))
      poly = sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(coords)), "x")))
      grid <- raster::raster(raster::extent(poly), nrows = ceiling(sqrt(nbins)), ncols= ceiling(sqrt(nbins)))
      sp::proj4string(grid)<-sp::proj4string(poly)
      
      cellLocationOnGrid = raster::cellFromXY(grid, currXY)
      tab <- table(cellLocationOnGrid)
      grid[as.numeric(names(tab))] <- tab
    }
    list(grid = grid, clusterIndexes = clusterIndexes, cellLocationOnGrid = cellLocationOnGrid, nFullbins = length(tab), maxBinSize = max(tab))
  })
  
  numOfRuns = ceiling(minSelection*max(unlist(lapply(initialGrids,function(clusterData){ clusterData$nFullbins * clusterData$maxBinSize  }))) / k)
  
  meanDistMatrix = rep(1,length(refNames))
  chosenCellList = lapply(1:numOfRuns, function(runNum){
    chosenCells = as.numeric(unlist(lapply(unique(refNames),function(currCluster){
      initialGrid = initialGrids[[which(unique(refNames)==currCluster)]]
      clusterIndexes = initialGrid$clusterIndexes
      grid = initialGrid$grid
      cellLocationOnGrid = initialGrid$cellLocationOnGrid
      kToUse = k
      if(k>length(which(!is.na(grid[])))){
        kToUse = length(which(!is.na(grid[])))
      }
      gridCellsToUse = sample(which(!is.na(grid[])),kToUse,replace = F)
      chosenCellsForCluster = clusterIndexes[unlist(lapply(gridCellsToUse, function(currCell){
        chosenCell = which(cellLocationOnGrid==currCell)
        if(length(chosenCell)>1){
          chosenCell = sample(chosenCell,1,prob = meanDistMatrix[clusterIndexes[chosenCell]])
        }
        chosenCell
      }))]
      chosenCellsForCluster
    })))
    cellsToReduce = chosenCells[which(meanDistMatrix[chosenCells]>minValueToReduceTo)]
    meanDistMatrix[cellsToReduce] <<- meanDistMatrix[cellsToReduce]/10
    chosenCells
  })
  
  chosenNeigList = lapply(1:length(refNames),function(cellIndex){
    selectedCellType = refNames[cellIndex]
    selectedCellIndexes = which(refNames == selectedCellType)
    cellXY = XY[cellIndex,]
    cellDist = fields::rdist(t(as.matrix(cellXY)),XY[selectedCellIndexes,])
    chosenRepeats = order(as.numeric(cellDist),decreasing = F)[1:neighborhoodSize]
    chosenRepeats = chosenRepeats[!is.na(chosenRepeats)]
    selectedCellIndexes[chosenRepeats]
  })
  list(chosenCellList = chosenCellList, chosenNeigList = chosenNeigList, numOfRuns = numOfRuns)
}
