# Finding Score for each Exon
findScore <- function(x)
{
    message("Performing Statistical testing of Gene")
    # Finding Score
    x <- cbind(x, 'Y' = apply(x, 1, function(x)(-log(as.numeric(x['P.Value'])))
                         *sign(as.numeric(x['FC']))))
    x <- cbind(x, 'Score' = apply(x, 1, function(x)
        (exp(-abs(as.numeric(x['Y']))))))
    x <- x[order(x$ID),]
    id <- x$ID
    id <- strsplit(id, ":", fixed=TRUE)
    x <- cbind(ID=as.character(lapply(id, FUN=function(x){paste(x[[1]][1])})),
               x)
    message('Aggregating Results')
    # Aggregating Results
    mExon <- cbind(aggregate(x$Score, list(x$ID), median),
                 aggregate(as.numeric(x$Score), list(x$ID), length)[2],
                 aggregate(as.numeric(x$FC), list(x$ID), median)[2],
                 aggregate(as.numeric(x$logFC), list(x$ID), median)[2])
    colnames(mExon) <- c('Gene', 'Median', 'ExonCount', 'FC', 'logFC')
    return(mExon)
}


