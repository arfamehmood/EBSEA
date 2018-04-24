# function to perform Exon Based Strategy
EBSEA <- function(data, group, paired = FALSE, plot = FALSE) {
  message("Checking Parameters")
  # Checking the parameters
  # Checking countData
  if(is.null(data)) {
    stop('No count data found..')
  }
  # Checking groups
  f <- factor(group)
  if(length(levels(f)) > 2) {
    stop('EBS performs pairwise comparison')
  }
  if(length(levels(f)) == 1) {
    stop('Please provide two groups')
  }
  if(length(levels(f)) == 0) {
    stop('No group information provided')
  }
  # Normalization

  norm.data <- normalizeData(data, group)

  # Designing Matrix
  if(paired == FALSE) {
    design <- model.matrix(~0 + f)
    colnames(design) <- levels(f)
    }
  if(paired == TRUE) {
    if(ncol(data) %% 2 != 0) {
      stop('Uneven number of columns found ')
    }
    if(length(which(f == levels(f)[1])) != length(which(f == levels(f)[2]))) {
      stop('Uneven samples in the data ....')
    }
    effects <- as.factor(rep(paste0('P', c(1:(ncol(data)/2))), 2))
    design <- model.matrix(~0 + f + effects)
    colnames(design) <- c(levels(f), levels(effects)[-1])
  }

  #Contrast matrix
  contrast <- paste(levels(f)[1], '-', levels(f)[2], sep = '')
  contrast.matrix <- makeContrasts(contrasts = contrast, levels = design)

  message("Performing Statistical testing of Exons")

  # Statistical testing
  voom.data <- voom(norm.data, design, plot = FALSE)
  fit <- lmFit(voom.data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit3 <- eBayes(fit2)
  exon.table <- topTable(fit3,
                         coef = 1,
                         adjust.method = "fdr",
                         number = nrow(norm.data$counts),
                         sort.by = "P")
  exon.table$ID <- rownames(exon.table)
  exon.table <- exon.table[, c("ID", "logFC", "P.Value", "adj.P.Val", "AveExpr")]
  #exon.table <- data.frame(cbind("FC" = logratio2foldchange(exon.table[["logFC"]]),
   #                              exon.table))
  # Generating Gene Statistics
  gene.table <- generateNullDistribution(findScore(exon.table))
  exon.table <- exon.table[, c(1, 5, 2, 3, 4)]
  colnames(exon.table) <- c('GeneExon', 'AveExpr', 'logFC',
                             'P.Value', 'FDR')
  row.names(exon.table) <- NULL

  result <- list("ExonTable" = exon.table,
                 "GeneTable" = gene.table,
                 "ExonCounts" = data,
                 "NormData" = norm.data$samples)

  tryCatch(if(plot == TRUE) {
    with(result$GeneTable, plot(logFC, -log2(P.Value),
        pch = 20, main = "Volcano plot"))
    abline(v = c(1, -1), lty = 3, lwd = 1)
  }, error = function(e){print(e)})
  class(result) <- "EBSEA"
  rm(exon.table)
  rm(gene.table)
  rm(norm.data)
  gc()
  message('Done')
  return(result)
}
