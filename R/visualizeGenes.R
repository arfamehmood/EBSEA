visualizeGenes <- function(gene, ebsea.out) {
    if(is.null(ebsea.out$ExonCounts)) {
      stop('No count data provided....')
    }
  group <- as.factor(ebsea.out$NormData$group)
  exon.table <- ebsea.out$ExonTable
  row.names(exon.table) <- exon.table$GeneExon
  # Fetching Information
  gene.exons <- ebsea.out$ExonTable[grep(gene, exon.table$GeneExon, fixed = TRUE), ]
  if(nrow(gene.exons) == 0) {
    stop('The gene provided is not in the list')
  }
  message('Obtaining Data')
  counts <- ebsea.out$ExonCounts[grep(gene, row.names(ebsea.out$ExonCounts), fixed = TRUE), ]
  counts <- t(apply(counts, 1, function(x)
          x * ebsea.out$NormData$norm.factors))
  gene.info <- ebsea.out$GeneTable[grep(gene, ebsea.out$GeneTable$Gene, fixed = TRUE), ]
  gene.exons <- gene.exons[order(gene.exons$GeneExon), ]

  # taking the mean and standard error of the rows of both samples group

  if(nrow(counts) > 1){
  mean.count <- data.frame('Group1' = apply(counts[, which(group %in% levels(group)[1])], 1, mean),
                           'Group2' = apply(counts[, which(group %in% levels(group)[2])], 1, mean))
  row.names(mean.count) <- c(1 : nrow(mean.count))
  se.count <- data.frame('Group1' = apply(counts[, group %in% levels(group)[1]],
                                    1, function(x)sd(x)/sqrt(length(x))),
                         'Group2' = apply(counts[, group %in% levels(group)[2]],
                                         1, function(x)sd(x)/sqrt(length(x))))
  } else{
      mean.count <- data.frame('Group1' = mean(counts[, which(group %in% levels(group)[1])]),
                               'Group2' = mean(counts[, which(group %in% levels(group)[2])]))

      se.count <- data.frame('Group1' = sd(counts[, group %in% levels(group)[1]])/sqrt(length(counts[, group %in% levels(group)[1]])),
                             'Group2' = sd(counts[, group %in% levels(group)[2]])/sqrt(length(counts[, group %in% levels(group)[2]])))
  }

  # Creating the plot
  message('Creating Plot')

  par(fig = c(0, 1, 0, 0.7), mar = c(4, 4, 2, 1), mgp = c(3, 1, 0))
  #plotting barplot of countData

  barx <- barplot(as.matrix(t(mean.count)), xlab = "Exon", ylab = "mean(Counts)",
                  col = c("gray", "white"), beside = TRUE,
                  ylim = c(0, max(mean.count) + max(se.count)), axis.lty = 1)

  toplotSE <- t(mean.count) + t(se.count) > (0.1 * sd(as.numeric(t(mean.count))))
  mask <- which(toplotSE)
  arrows(barx[mask], t(mean.count)[mask] + t(se.count)[mask], barx[mask],
  t(mean.count)[mask] - t(se.count)[mask], angle = 90,
           code = 3, length = 0.03)
  par(fig = c(0, 1, 0, 0.7))
  add_legend("bottomright", legend = levels(group), fill = c("grey", "white"),
         horiz = TRUE, cex = 0.9, box.lty = 0)
  par(fig = c(0, 1, 0.7, 0.98), new = TRUE, mar = c(0.5, 4, 2, 1))
  FDR <- c()
  for( i in 1 : nrow(gene.exons)) {
      if(gene.exons$P.Value[i] < 0.01) {
          FDR <- c(FDR, '**')
      } else if(gene.exons$P.Value[i] < 0.05) {
          FDR <- c(FDR, '*')
      } else  {
          FDR <- c(FDR, ' ')
      }
  }

  if(max(gene.exons$logFC) <= 0){
        max.value <- 0
  } else if(max(gene.exons$logFC) > 0 & max(gene.exons$logFC) <= 1){
      max.value <- max(gene.exons$logFC) + 0.1
  } else {
  max.value <- max(gene.exons$logFC) + 0.5
  }
  if(min(gene.exons$logFC) <= 0){
      min.value <- min(gene.exons$logFC) - 0.5
  } else if(min(gene.exons$logFC) < 0 & min(gene.exons$logFC) >= -1) {
      min.value <- min(gene.exons$logFC) - 0.1
  } else{
      min.value <- 0
  }
  b2 <- barplot(gene.exons$logFC, col = ifelse(gene.exons$logFC > 0, "red", "blue"),
                 ylab = 'Log FC',
                 ylim = c(min.value, max.value))

   text(b2, y = ifelse(gene.exons$logFC < 0,
                           gene.exons$logFC - 0.1,
                           gene.exons$logFC + 0.1),
            labels = FDR)
   gene <- strsplit(row.names(counts)[1], ':')[[1]][1]
   par(fig = c(0, 1, 0.98, 1), new = TRUE)
   mtext(paste0(gene, ' (FDR:', round(gene.info$FDR, 5), ')'))
   par(fig = c(0, 1, 0, 1))
}


add_legend <- function(...) {
  opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),
                mar = c(0, 0, 0, 0), new = TRUE)
    on.exit(par(opar))
    plot(0, 0, type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n')
    legend(...)
}
