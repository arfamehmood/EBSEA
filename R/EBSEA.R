# function to perform Exon Based Strategy
EBSEA <- function(countData, group, paired=FALSE, effects=NULL, plot=FALSE)
{
    message("Checking Parameters")

    # Checking the parameters
    # Checking countData
    if(is.null(countData)){stop('No count data found..')}
    # Checking groups
    f <- factor(group)
    if(length(levels(f)) > 2){stop('EBS performs pairwise comparison')}
    if(length(levels(f)) == 1){stop('Please provide two groups')}
    if(length(levels(f)) == 0){stop('No group information provided')}
    # Normalization
    countData <- normalizeData(countData, group)

    # Designing Matrix
    if(paired == FALSE)
    {
        design <- model.matrix(~0 + f)
        colnames(design) <- levels(f)
    }
    if(paired == TRUE)
    {
        effects <- as.factor(effects)
        design <- model.matrix(~0 + f + effects)
        colnames(design) <- c(levels(f), levels(effects)[-1])
    }

    # Contrast matrix
    contrast <- paste(levels(f)[1], '-', levels(f)[2], sep='')
    contrast.matrix <- makeContrasts(contrasts=contrast, levels=design)

    message("Performing Statistical testing of Exons")

    # Statistical testing
    v <- voom(countData, design, plot=FALSE)
    fit <- lmFit(v, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit3 <- eBayes(fit2)
    ExonTable <- topTable(fit3, coef=1, adjust.method="fdr",
                          number=nrow(countData$counts), sort.by="P")

    ExonTable$ID <- rownames(ExonTable)
    ExonTable <- ExonTable[, c("ID", "logFC", "P.Value", "adj.P.Val", "AveExpr")]

    ExonTable <- data.frame(cbind("FC" = logratio2foldchange(ExonTable[["logFC"]]),
                                  ExonTable))

    # Generating Gene Statistics
    GeneTable <- generateNullDistribution(findScore(ExonTable))
    ExonTable <- ExonTable[, c(2, 6, 1, 3, 4, 5)]
    colnames(ExonTable) <- c('GeneExon', 'AveExpr', 'FC', 'logFC',
                             'P.Value', 'FDR')
    row.names(ExonTable) <- NULL


    result <- list("ExonTable" = ExonTable,
                 "GeneTable" = GeneTable)

    if(plot == TRUE)
    {
        with(result$GeneTable, plot(logFC, -log2(P.Value), pch=20,
                                    main="Volcano plot"))
        abline(v=c(1, -1), lty=3, lwd=1)
    }

    rm(ExonTable)
    gc()
    message('Done')
    return(result)

}
