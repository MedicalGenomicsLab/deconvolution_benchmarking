library(quadprog)
library(reshape)
library(e1071)
library(Seurat)
library(ROCR)
library(varhandle)
library(MAST)


# perform DE analysis using Seurat
DEAnalysis <- function(scdata, i, id, path) {
    exprObj <- CreateSeuratObject(raw.data = as.data.frame(scdata), project = "DE")
    exprObj2 <- SetIdent(exprObj, ident.use = as.vector(id))
    print("Calculating differentially expressed genes:")
    
    de_group <- FindMarkers(
        object = exprObj2, ident.1 = i, ident.2 = NULL,
        only.pos = TRUE, test.use = "bimod"
    )
    save(de_group, file = paste(path, "/de_", i, ".RData", sep = ""))
}

# build signature matrix using genes identified by DEAnalysis()
buildSignatureMatrixUsingSeurat <- function(scdata, id, path, diff.cutoff = 0.5, pval.cutoff = 0.01) {
    # # perform differential expression analysis
    # DEAnalysis(scdata, id, path)

    numberofGenes <- c()
    for (i in unique(id)) {
        load(file = paste(path, "/de_", i, ".RData", sep = ""))
        DEGenes <- rownames(de_group)[intersect(which(de_group$p_val_adj < pval.cutoff), which(de_group$avg_logFC > diff.cutoff))]
        nonMir <- grep("MIR|Mir", DEGenes, invert = T)
        assign(paste("cluster_lrTest.table.", i, sep = ""), de_group[which(rownames(de_group) %in% DEGenes[nonMir]), ])
        numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
    }

    # need to reduce number of genes
    # for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
    # choose matrix with lowest condition number
    conditionNumbers <- c()
    for (G in 50:200) {
        Genes <- c()
        j <- 1
        for (i in unique(id)) {
            if (numberofGenes[j] > 0) {
                temp <- paste("cluster_lrTest.table.", i, sep = "")
                temp <- as.name(temp)
                temp <- eval(parse(text = temp))
                temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
                Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
            }
            j <- j + 1
        }
        Genes <- unique(Genes)
        # make signature matrix
        ExprSubset <- scdata[Genes, ]
        Sig <- NULL
        for (i in unique(id)) {
            Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
        }
        colnames(Sig) <- unique(id)
        conditionNumbers <- c(conditionNumbers, kappa(Sig))
    }
    G <- which.min(conditionNumbers) + min(49, numberofGenes - 1) # G is optimal gene number
    #
    Genes <- c()
    j <- 1
    for (i in unique(id)) {
        if (numberofGenes[j] > 0) {
            temp <- paste("cluster_lrTest.table.", i, sep = "")
            temp <- as.name(temp)
            temp <- eval(parse(text = temp))
            temp <- temp[order(temp$p_val_adj, decreasing = TRUE), ]
            Genes <- c(Genes, (rownames(temp)[1:min(G, numberofGenes[j])]))
        }
        j <- j + 1
    }
    Genes <- unique(Genes)
    ExprSubset <- scdata[Genes, ]
    Sig <- NULL
    for (i in unique(id)) {
        Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
    }
    colnames(Sig) <- unique(id)
    save(Sig, file = paste(path, "/Sig.RData", sep = ""))
    return(Sig)
}

## alternative differential expression method using MAST

# functions for DE

Mean.in.log2space <- function(x, pseudo.count) {
    return(log2(mean(2^(x) - pseudo.count) + pseudo.count))
}

stat.log2 <- function(data.m, group.v, pseudo.count) {
    # data.m=data.used.log2
    log2.mean.r <- aggregate(t(data.m), list(as.character(group.v)), function(x) Mean.in.log2space(x, pseudo.count))
    log2.mean.r <- t(log2.mean.r)
    colnames(log2.mean.r) <- paste("mean.group", log2.mean.r[1, ], sep = "")
    log2.mean.r <- log2.mean.r[-1, ]
    log2.mean.r <- as.data.frame(log2.mean.r)
    log2.mean.r <- varhandle::unfactor(log2.mean.r) # from varhandle
    log2.mean.r[, 1] <- as.numeric(log2.mean.r[, 1])
    log2.mean.r[, 2] <- as.numeric(log2.mean.r[, 2])
    log2_foldchange <- log2.mean.r$mean.group1 - log2.mean.r$mean.group0
    results <- data.frame(cbind(log2.mean.r$mean.group0, log2.mean.r$mean.group1, log2_foldchange))
    colnames(results) <- c("log2.mean.group0", "log2.mean.group1", "log2_fc")
    rownames(results) <- rownames(log2.mean.r)
    return(results)
}

v.auc <- function(data.v, group.v) {
    prediction.use <- prediction(data.v, group.v, 0:1)
    perf.use <- performance(prediction.use, "auc")
    auc.use <- round(perf.use@y.values[[1]], 3)
    return(auc.use)
}
m.auc <- function(data.m, group.v) {
    AUC <- apply(data.m, 1, function(x) v.auc(x, group.v))
    AUC[is.na(AUC)] <- 0.5
    return(AUC)
}

# perform DE analysis using MAST
DEAnalysisMAST <- function(scdata, i, id, path) {
    pseudo.count <- 0.1
    data.used.log2 <- log2(scdata + pseudo.count)
    colnames(data.used.log2) <- make.unique(colnames(data.used.log2))
    diff.cutoff <- 0.5

    cells.symbol.list2 <- colnames(data.used.log2)[which(id == i)]
    cells.coord.list2 <- match(cells.symbol.list2, colnames(data.used.log2))
    cells.symbol.list1 <- colnames(data.used.log2)[which(id != i)]
    cells.coord.list1 <- match(cells.symbol.list1, colnames(data.used.log2))
    data.used.log2.ordered <- cbind(
        data.used.log2[, cells.coord.list1], data.used.log2[, cells.coord.list2]
    )
    group.v <- c(
        rep(0, length(cells.coord.list1)), rep(1, length(cells.coord.list2))
    )
    # ouput
    log2.stat.result <- stat.log2(data.used.log2.ordered, group.v, pseudo.count)
    Auc <- m.auc(data.used.log2.ordered, group.v)
    bigtable <- data.frame(cbind(log2.stat.result, Auc))

    DE <- bigtable[bigtable$log2_fc > diff.cutoff, ]
    dim(DE)
    if (dim(DE)[1] > 1) {
        data.1 <- data.used.log2[, cells.coord.list1]
        data.2 <- data.used.log2[, cells.coord.list2]
        genes.list <- rownames(DE)
        log2fold_change <- cbind(genes.list, DE$log2_fc)
        colnames(log2fold_change) <- c("gene.name", "log2fold_change")
        counts <- as.data.frame(
            cbind(data.1[genes.list, ], data.2[genes.list, ])
        )
        groups <- c(
            rep(
                "Cluster_Other",
                length(cells.coord.list1)
            ),
            rep(i, length(cells.coord.list2))
        )
        groups <- as.character(groups)
        data_for_MIST <- as.data.frame(
            cbind(
                rep(rownames(counts), dim(counts)[2]),
                melt(counts),
                rep(groups, each = dim(counts)[1]),
                rep(1, dim(counts)[1] * dim(counts)[2])
            )
        )
        colnames(data_for_MIST) <- c(
            "Gene", "Subject.ID", "Et", "Population", "Number.of.Cells"
        )
        vbeta <- data_for_MIST
        vbeta.fa <- FromFlatDF(vbeta,
            idvars = c("Subject.ID"),
            primerid = "Gene", measurement = "Et", ncells = "Number.of.Cells",
            geneid = "Gene", cellvars = c("Number.of.Cells", "Population"),
            phenovars = c("Population"), id = "vbeta all"
        )
        vbeta.1 <- subset(vbeta.fa, Number.of.Cells == 1)
        # .3 MAST
        head(colData(vbeta.1))
        zlm.output <- zlm(
            ~Population, vbeta.1,
            method = "bayesglm", ebayes = TRUE
        )
        show(zlm.output)
        coefAndCI <- summary(zlm.output, logFC = TRUE)
        zlm.lr <- lrTest(zlm.output, "Population")
        zlm.lr_pvalue <- melt(zlm.lr[, , "Pr(>Chisq)"])
        zlm.lr_pvalue <- zlm.lr_pvalue[
            which(zlm.lr_pvalue$test.type == "hurdle"),
        ]



        lrTest.table <- merge(
            zlm.lr_pvalue, DE,
            by.x = "primerid", by.y = "row.names"
        )
        colnames(lrTest.table) <- c(
            "Gene",
            "test.type",
            "p_value",
            paste("log2.mean.", "Cluster_Other", sep = ""),
            paste("log2.mean.", i, sep = ""),
            "log2fold_change",
            "Auc"
        )
        cluster_lrTest.table <- lrTest.table[rev(order(lrTest.table$Auc)), ]

        # . 4 save results
        write.csv(
            cluster_lrTest.table,
            file = paste(path, "/", i, "_lrTest.csv", sep = "")
        )
        save(
            cluster_lrTest.table,
            file = paste(path, "/", i, "_MIST.RData", sep = "")
        )
    }
}

# build signature matrix using genes identified by DEAnalysisMAST()
buildSignatureMatrixMAST <- function(scdata, id, path, diff.cutoff = 0.5, pval.cutoff = 0.01) {
    # compute differentially expressed genes for each cell type
    # DEAnalysisMAST(scdata, id, path)

    # for each cell type, choose genes in which FDR adjusted p-value is less than 0.01 and the estimated fold-change
    # is greater than 0.5
    numberofGenes <- c()
    for (i in unique(id)) {
        if (file.exists(paste(path, "/", i, "_MIST.RData", sep = ""))) {
            load(file = paste(path, "/", i, "_MIST.RData", sep = ""))
            pvalue_adjusted <- p.adjust(cluster_lrTest.table[, 3], method = "fdr", n = length(cluster_lrTest.table[, 3]))
            cluster_lrTest.table <- cbind(cluster_lrTest.table, pvalue_adjusted)
            DEGenes <- cluster_lrTest.table$Gene[intersect(which(pvalue_adjusted < pval.cutoff), which(cluster_lrTest.table$log2fold_change > diff.cutoff))]
            nonMir <- grep("MIR|Mir", DEGenes, invert = T) # because Mir gene is usually not accurate
            assign(paste("cluster_lrTest.table.", i, sep = ""), cluster_lrTest.table[which(cluster_lrTest.table$Gene %in% DEGenes[nonMir]), ])
            numberofGenes <- c(numberofGenes, length(DEGenes[nonMir]))
        }
    }

    # need to reduce number of genes
    # for each subset, order significant genes by decreasing fold change, choose between 50 and 200 genes
    # for each, iterate and choose matrix with lowest condition number
    conditionNumbers <- c()
    for (G in 50:200) {
        Genes <- c()
        j <- 1
        for (i in unique(id)) {
            if (numberofGenes[j] > 0) {
                temp <- paste("cluster_lrTest.table.", i, sep = "")
                temp <- as.name(temp)
                temp <- eval(parse(text = temp))
                temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
                Genes <- c(Genes, varhandle::unfactor(temp$Gene[1:min(G, numberofGenes[j])]))
            }
            j <- j + 1
        }
        Genes <- unique(Genes)
        # make signature matrix
        ExprSubset <- scdata[Genes, ]
        Sig <- NULL
        for (i in unique(id)) {
            Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
        }
        colnames(Sig) <- unique(id)
        conditionNumbers <- c(conditionNumbers, kappa(Sig))
    }
    G <- which.min(conditionNumbers) + min(49, numberofGenes - 1)
    Genes <- c()
    j <- 1
    for (i in unique(id)) {
        if (numberofGenes[j] > 0) {
            temp <- paste("cluster_lrTest.table.", i, sep = "")
            temp <- as.name(temp)
            temp <- eval(parse(text = temp))
            temp <- temp[order(temp$log2fold_change, decreasing = TRUE), ]
            Genes <- c(Genes, varhandle::unfactor(temp$Gene[1:min(G, numberofGenes[j])]))
        }
        j <- j + 1
    }
    Genes <- unique(Genes)
    ExprSubset <- scdata[Genes, ]
    Sig <- NULL
    for (i in unique(id)) {
        Sig <- cbind(Sig, (apply(ExprSubset, 1, function(y) mean(y[which(id == i)]))))
    }
    colnames(Sig) <- unique(id)
    save(Sig, file = paste(path, "/Sig.RData", sep = ""))
    return(Sig)
}
