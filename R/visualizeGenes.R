#' Visualize gene
#'
#' Produces a visualization summarizing the normalized read count
#' in each sample group and expression difference across the
#' expressed exons.Top panel contains the log2 fold-change for
#' each expressed exon. Asterisk denotes the significance
#' level (*: < 0.05, **: < 0.01). Bottom panel shows the averaged normalized
#' read count for each sample group. The title of the figure shows
#' the gene name and the adjusted gene-level p-value (fdr)
#' @param gene Gene name matching the input data.
#' @param ebsea.out Plots the mean count and fold-change
#' the exons of the specified gene.
#' @return Plots the mean count and fold-change across the
#' exons of the specified gene.
#' @examples data(exonCounts)
#' @examples group <- data.frame('group' = as.factor(c('G1', 'G1', 'G1', 'G2', 'G2', 'G2', 'G2')))
#' @examples row.names(group) <- colnames(exonCounts)
#' @examples design <- ~group
#' @examples ebsea.out <- EBSEA(exonCounts, group, design)
#' @examples visualizeGenes('FBgn0000017', ebsea.out)
#' @importFrom graphics abline mtext par text arrows barplot legend
#' @importFrom stats sd
#' @export visualizeGenes



visualizeGenes <- function(gene, ebsea.out) {
  if(is.null(ebsea.out$RawCounts)) {
    stop('No count data provided....')
  }
  group <- as.factor(ebsea.out$Group)
  exon.table <- ebsea.out$ExonTable
  exon.table$GeneExon <- row.names(exon.table)
  # Fetching Information
  gene.exons <- exon.table[grep(gene, exon.table$GeneExon, fixed = TRUE), , drop = FALSE]
  if(nrow(gene.exons) == 0) {
    stop('The gene provided is not in the list')
  }
  message('Obtaining Data')
  counts <- ebsea.out$NormCounts[grep(gene, row.names(ebsea.out$NormCounts),
                                      fixed = TRUE), , drop = FALSE]
  gene.info <- ebsea.out$GeneTable[grep(gene, ebsea.out$GeneTable$Gene, fixed = TRUE),
                                   , drop = FALSE]
  gene.exons <- gene.exons[order(gene.exons$GeneExon), ]
  
  # taking the mean and standard error of the rows of both samples
  mean.count <- data.frame('Group1' = apply(counts[, which(group %in% levels(group)[1]), drop = FALSE], 1, mean),
                           'Group2' = apply(counts[, which(group %in% levels(group)[2]), drop = FALSE], 1, mean))
  se.count <- data.frame('Group1' = apply(counts[, group %in% levels(group)[1], drop = FALSE],
                                          1, function(x)sd(x)/sqrt(length(x))),
                         'Group2' = apply(counts[, group %in% levels(group)[2], drop = FALSE],
                                          1, function(x)sd(x)/sqrt(length(x))))
  # checking if the p-values have NA values. If NA is found it is converted to 1
  gene.exons$pvalue[is.na(gene.exons$pvalue)] <- 1
  
  # Creating the plot
  message('Creating Plot')
  #plot.new()
  #frame()
  par(fig=c(0, 1, 0, 0.7), mar = c(4, 4, 0.5, 1), mgp = c(3, 1, 0))
  #plotting barplot of countData
  barx <- barplot(as.matrix(t(mean.count)), xlab = "Exons", ylab = "mean(Counts)",
                  col = c("gray", "white"),  beside = TRUE,
                  ylim = c(0, max(mean.count) + max(se.count)),  names.arg =seq(1, nrow(mean.count)))
  toplotSE <- t(mean.count) + t(se.count) > (0.1 * sd(as.numeric(t(mean.count))))
  mask <- which(toplotSE)
  arrows(barx[mask], t(mean.count)[mask] + t(se.count)[mask], barx[mask],
         t(mean.count)[mask] - t(se.count)[mask], angle = 90,
         code = 3, length = 0.03)
  par(fig = c(0, 1, 0, 0.7))
  add_legend("bottomright", legend = levels(group), fill = c("grey", "white"),
             horiz = TRUE, cex = 0.9, box.lty = 0)
  par(fig = c(0, 1, 0.725, 1.0), new = TRUE, mar = c(0.5, 4, 2, 1))
  
  # setting the arrows for the bars based on their p-values
  FDR <- c()
  for( i in seq(1, nrow(gene.exons))) {
    if(gene.exons$pvalue[i] <= 0.001) {
      FDR <- c(FDR, '***')
    } else if(gene.exons$pvalue[i] <= 0.01 & gene.exons$pvalue[i] > 0.001) {
      FDR <- c(FDR, '**')
    } else if(gene.exons$pvalue[i] <= 0.05 & gene.exons$pvalue[i] > 0.01) {
      FDR <- c(FDR, '*')
    }else  {
      FDR <- c(FDR, '')
    }
  }
  
  # finding the upper and lower limit of the log2FC plot
  sl <- ifelse(min(gene.exons$log2FoldChange) < 0, floor(min(gene.exons$log2FoldChange)), 0)
  sl <- ifelse(min(gene.exons$log2FoldChange) > 0, 0, sl )
  el <- ifelse(max(gene.exons$log2FoldChange) < 0, floor(max(gene.exons$log2FoldChange)), 0 )
  el <- ifelse(max(gene.exons$log2FoldChange) > 0, ceiling(max(gene.exons$log2FoldChange)), 0)
  if(round(min(gene.exons$log2FoldChange)) == sl){
    sl <- sl - 0.5
  } else if(round(max(gene.exons$log2FoldChange)) == el){
    el <- el + 0.5
  }
  b2 <- barplot(gene.exons$log2FoldChange, col = ifelse(gene.exons$log2FoldChange > 0, "red", "blue"),
                ylab = 'Log2 FC',
                ylim = c(sl, el))
  
  
  text(b2, y =  ifelse(gene.exons$log2FoldChange < 0,
                       gene.exons$log2FoldChange + (sl + abs(gene.exons$log2FoldChange))/3,
                       gene.exons$log2FoldChange + (el - gene.exons$log2FoldChange)/3),
       labels = FDR, srt = 90)
  gene <- strsplit(row.names(counts)[1], ':')[[1]][1]
  # mtext(paste0(gene, ' (FDR:', round(gene.info$padj, 3), ')'))
  mtext(paste0(gene, ' (FDR:', round(gene.info$padj, 3), ')'), line  = 0.5)
}


add_legend <- function(...) {
  opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
              mar = c(0, 0, 0, 0), new = TRUE)
  on.exit(par(opar))
  plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend(...)
}
