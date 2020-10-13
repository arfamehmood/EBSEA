#' Filter Count Data
#'
#' Filtering of exons based on their expression levels
#' @param x A numeric dataframe of exon counts across the samples.
#' Exon number in format GeneName:Exonnumber should be indicated
#' in the row name and sample names as column names.
#' @param  mean Exons with average count value across the
#' dataset less than mean are filtered out. Default: 1
#' @param exonCount After filtering the individual exons,
#' only genes with at least the given number of exons
#' remaining will be retained. Default: 1
#' @return A dataframe of filtered counts of exons
#' @examples
#' data(exonCounts)
#' res <- filterCounts(exonCounts)
#' @export filterCounts
filterCounts <- function(x, mean = 1, exonCount = 1) {
  if(ncol(x) == 0 | nrow(x) == 0) {
    stop('No count data found')
  } else {
    message(paste0('Filter exons with average count less than ',
                   mean, ' across the dataset and keep genes with
                     at least ', exonCount  ,' remaining exon'))
    x <- x[rowMeans(as.matrix(x)) > mean, ]
    ## getting the exon count of each gene
    x$genes <- gsub(':.*', '', row.names(x))
    genes <- x$genes[x$genes >= exonCount]
    x <- x[x$genes %in% genes, ]
    
    return(x[, seq(1, (ncol(x)-1))])
  }
}
