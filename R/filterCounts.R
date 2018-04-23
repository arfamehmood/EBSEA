filterCounts <- function(x, mean=1) {
  if(ncol(x) == 0 | nrow(x) == 0) {
    stop('No count data found')
  } else {
    message('Filtering genes having normalized mean less than 1')
    normalize.count <- normalizeData(x)
    counts <- t(apply(x, 1, function(y)
    y*normalize.count$samples$norm.factors))
    counts <- cbind(counts, 'rowMeans' = rowSums(counts)/ncol(counts))
    counts <- cbind(x,counts)
    counts <- counts[counts[, ncol(counts)] > mean, ]
    return(x[,1:ncol(x)-1])
    }
}
