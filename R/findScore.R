# Finding Score for each Exon
findScore <- function(x) {
  message("Performing Statistical testing of Gene")
  # Finding Score
  x <- cbind(x, 'Y' = apply(x, 1, function(x)(-log(as.numeric(x['P.Value'])))
                      *sign(as.numeric(x['logFC']))))
  x <- x[order(x$ID),]
  id <- x$ID
  id <- strsplit(id, ":", fixed = TRUE)
  x <- cbind('Gene' = as.character(lapply(id, FUN = function(x){paste(x[[1]][1])})),
               x)
  median.score <-  aggregate(x$Y, list(x$Gene), median)
  colnames(median.score)[1] <- 'Gene'
  x <- join(median.score, x , by = 'Gene', type = 'right')

  x <- cbind(x, 'Score' = apply(x, 1, function(x)
     (exp(-abs(as.numeric(x['x']))))))
  message('Aggregating Results')
  # Aggregating Results
  median.exon <- cbind(aggregate(x$Score, list(x$Gene), median),
                 aggregate(as.numeric(x$Score), list(x$Gene), length)[2],
                # aggregate(as.numeric(x$FC), list(x$Gene), median)[2],
                 aggregate(as.numeric(x$logFC), list(x$Gene), median)[2])
  colnames(median.exon) <- c('Gene', 'Median(signed_P.Value)', 'ExonCount', 'logFC')
  return(median.exon)
}
