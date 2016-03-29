filterCounts <- function(x, noOfSamples)
{

    if(ncol(x) == 0 | nrow(x) == 0)
    {
        stop('No count data found')
    }
    else
    {
        message(paste('Filtering genes having a cpm of more than 1 in atleast',
                      noOfSamples, 'percent of the sample.'))
        n <- round((ncol(x) * (100-noOfSamples))/100, digits=0)
        isexpr <- rowSums(cpm(x) > 1) >= n
        x <- x[isexpr, ]
        return(x)
    }
}
