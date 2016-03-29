visualizeGenes <- function(gene, group, countData, result)
{
    if(is.null(countData)){stop('No count data provided....')}
    group <- as.factor(group)
    ExonTable <- result$ExonTable
    row.names(ExonTable) <- ExonTable[, 1]

    # Fetching Information
    geneExons <- result$ExonTable[grep(gene, row.names(ExonTable)), ]
    if(nrow(geneExons) == 0)
    {stop('The gene provided is not in the list')}
    message('Obtaining Data')
    norCount <- normalizeData(countData, group)
    countData <- t(apply(countData, 1, function(x)
        x*norCount$samples$norm.factors))
    counts <- countData[grep(gene, row.names(countData)), ]
    geneInfo <- result$GeneTable[grep(gene, result$GeneTable$Gene), ]
    geneExons <- geneExons[order(geneExons$GeneExon), ]

    # taking the mean and standard error of the rows of both samples
    meanCount <- data.frame(apply(counts[, group %in% levels(group)[1]], 1, mean),
                        apply(counts[, group %in% levels(group)[2]], 1, mean))
    row.names(meanCount) <- c(1:nrow(meanCount))
    seCount<-data.frame('Group1' = apply(counts[, group %in% levels(group)[1]],
                                    1, function(x)sd(x)/length(x)),
                        'Group2' = apply(counts[, group %in% levels(group)[2]],
                                         1, function(x)sd(x)/length(x)))

    # Creating the plot
    message('Creating Plot')

   # plot.new()
    #frame()

    par(fig=c(0, 1, 0, 0.7), mar=c(4, 4, 0.5, 1), mgp=c(3, 1, 0))

    #plotting barplot of countData
    barx <- barplot(as.matrix(t(meanCount)), xlab="Exon", ylab="mean(Counts)",
                  col=c("gray", "white"), axis.lty = 1, beside=TRUE,
                  ylim=c(0, max(meanCount)+max(seCount)))
    toplotSE <- t(meanCount) + t(seCount) > (0.1 * sd(as.numeric(t(meanCount))))
    mask <- which(toplotSE)
    arrows(barx[mask], t(meanCount)[mask] + t(seCount)[mask], barx[mask],
           t(meanCount)[mask] - t(seCount)[mask], angle=90,
           code=3, length=0.03)
   par(fig=c(0, 1, 0, 0.7))
   add_legend("bottomright", legend=levels(group), fill=c("grey", "white"),
            horiz=TRUE, cex=0.9, box.lty=0)
   par(fig=c(0, 1, 0.7, 1), new=TRUE, mar=c(0.5, 4, 2, 1))
   FDR <- c()
   for( i in 1:nrow(geneExons))
   {if(geneExons$FDR[i] < 0.005){FDR <- c(FDR, '**')}
       else if(geneExons$FDR[i] < 0.05 & geneExons$FDR[i] >= 0.005)
       {FDR <- c(FDR, '*')}
       else  {FDR <- c(FDR, ' ')}
   }

   b2 <- barplot(geneExons$FC, col=ifelse(geneExons$FC > 0, "red", "blue"),
                 ylab='Fold Change', ylim=c(min(geneExons$FC) - 1.5,
                                            max(geneExons$FC) + 1.5))
   text(b2, y=ifelse(geneExons$FC < 0, geneExons$FC - 1, geneExons$FC + 1),
        labels=FDR)
   mtext(paste(gene, '( FDR:', round(geneInfo$FDR, 3), ')') )
}


add_legend <- function(...) {
    opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0),
                mar=c(0, 0, 0, 0), new=TRUE)
    on.exit(par(opar))
    plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
    legend(...)
}
