#calculating median
medianP <- function(x) {
    out <- median(x, na.rm = TRUE)
    out <- exp(-abs(out))# This ensures that the output are p-value type scores
    #(median p-value).
    out
}

# Generating p-values/FDR for each gene
generateNullDistribution <- function(x) {
  set.seed(1234)
  exon.counts <- unique(x[, 'ExonCount'])
  x$P.Value <- 0
  N <- 100000
  # Generating NULL Distribution and calculating p-values for the genes
  message('Generating Distribution and generating P-values')
  for(i in 1:length(exon.counts)) {
    null.distribution <- NULL
      for(j in 1:exon.counts[i]) {
        null.distribution <- cbind(null.distribution,
                                   sample(c(-1, 1),
            N, replace = TRUE) * (-log(runif(N, min = 0, max = 1))))
        }
        M <- apply(null.distribution, 1, medianP)
        ind <- which(x$ExonCount == exon.counts[i])
        med <- x$Median[ind]
        temp <- sapply(med, function(x)(sum(x >= M) + 1)/(N + 1))
        x$P.Value[ind] <- as.numeric(temp)
    }
    # FDR Correction
    message('Correcting P-values')
    x <- cbind(x, 'FDR' = p.adjust(x$P.Value, "BH"))
    rm(null.distribution)
    return(x)
}
