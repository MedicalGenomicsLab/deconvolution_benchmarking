# load packages
library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)

# trim bulk and single-cell data to contain the same genes
trimData <- function(Signature, bulkData) {
    Genes <- intersect(rownames(Signature), names(bulkData))
    B <- bulkData[Genes]
    S <- Signature[Genes, ]
    return(list("sig" = S, "bulk" = B))
}


# solve using OLS, constrained such that cell type numbers>0
solveOLS <- function(S, B) {
    D <- t(S) %*% S
    d <- t(S) %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    solution <- solve.QP(D, d, A, bzero)$solution
    names(solution) <- colnames(S)
    print(round(solution / sum(solution), 5))
    return(solution / sum(solution))
}

# return cell number, not proportion
# do not print output
solveOLSInternal <- function(S, B) {
    D <- t(S) %*% S
    d <- t(S) %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    solution <- solve.QP(D, d, A, bzero)$solution
    names(solution) <- colnames(S)
    return(solution)
}

# solve using WLS with weights dampened by a certain dampening constant
solveDampenedWLS <- function(S, B) {
    # first solve OLS, use this solution to find a starting point for the weights
    solution <- solveOLSInternal(S, B)
    # now use dampened WLS, iterate weights until convergence
    iterations <- 0
    changes <- c()
    # find dampening constant for weights using cross-validation
    j <- findDampeningConstant(S, B, solution)
    change <- 1
    while (change > .01 & iterations < 1000) {
        newsolution <- solveDampenedWLSj(S, B, solution, j)
        # decrease step size for convergence
        solutionAverage <- rowMeans(cbind(newsolution, matrix(solution, nrow = length(solution), ncol = 4)))
        change <- norm(as.matrix(solutionAverage - solution))
        solution <- solutionAverage
        iterations <- iterations + 1
        changes <- c(changes, change)
    }
    print(round(solution / sum(solution), 5))
    return(solution / sum(solution))
}

# solve WLS given a dampening constant
solveDampenedWLSj <- function(S, B, goldStandard, j) {
    multiplier <- 1 * 2^(j - 1)
    sol <- goldStandard
    ws <- as.vector((1 / (S %*% sol))^2)
    wsScaled <- ws / min(ws)
    wsDampened <- wsScaled
    wsDampened[which(wsScaled > multiplier)] <- multiplier
    W <- diag(wsDampened)
    D <- t(S) %*% W %*% S
    d <- t(S) %*% W %*% B
    A <- cbind(diag(dim(S)[2]))
    bzero <- c(rep(0, dim(S)[2]))
    sc <- norm(D, "2")
    solution <- solve.QP(D / sc, d / sc, A, bzero)$solution
    names(solution) <- colnames(S)
    return(solution)
}

# find a dampening constant for the weights using cross-validation
findDampeningConstant <- function(S, B, goldStandard) {
    solutionsSd <- NULL
    # goldStandard is used to define the weights
    sol <- goldStandard
    ws <- as.vector((1 / (S %*% sol))^2)
    wsScaled <- ws / min(ws)
    wsScaledMinusInf <- wsScaled
    # ignore infinite weights
    if (max(wsScaled) == "Inf") {
        wsScaledMinusInf <- wsScaled[-which(wsScaled == "Inf")]
    }
    # try multiple values of the dampening constant (multiplier)
    # for each, calculate the variance of the dampened weighted solution for a subset of genes
    for (j in 1:ceiling(log2(max(wsScaledMinusInf)))) {
        multiplier <- 1 * 2^(j - 1)
        wsDampened <- wsScaled
        wsDampened[which(wsScaled > multiplier)] <- multiplier
        solutions <- NULL
        seeds <- c(1:100)
        for (i in 1:100) {
            set.seed(seeds[i]) # make nondeterministic
            subset <- sample(length(ws), size = length(ws) * 0.5) # randomly select half of gene set
            # solve dampened weighted least squares for subset
            fit <- lm(B[subset] ~ -1 + S[subset, ], weights = wsDampened[subset])
            sol <- fit$coef * sum(goldStandard) / sum(fit$coef)
            solutions <- cbind(solutions, sol)
        }
        solutionsSd <- cbind(solutionsSd, apply(solutions, 1, sd))
    }
    # choose dampening constant that results in least cross-validation variance
    j <- which.min(colMeans(solutionsSd^2))
    return(j)
}

solveSVR <- function(S, B) {
    # scaling
    ub <- max(c(as.vector(S), B)) # upper bound
    lb <- min(c(as.vector(S), B)) # lower bound
    Bs <- ((B - lb) / ub) * 2 - 1
    Ss <- ((S - lb) / ub) * 2 - 1

    # perform SVR
    model <- svm(Ss, Bs, nu = 0.5, scale = TRUE, type = "nu-regression", kernel = "linear", cost = 1)
    coef <- t(model$coefs) %*% model$SV
    coef[which(coef < 0)] <- 0
    coef <- as.vector(coef)
    names(coef) <- colnames(S)
    print(round(coef / sum(coef), 5))
    return(coef / sum(coef))
}
