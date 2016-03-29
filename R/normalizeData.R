normalizeData <- function(x, group)
{
    message('Normalizing Data')
    y <- DGEList(counts=x, group=group)
    y <- calcNormFactors(y, method='TMM', na.rm=TRUE)
    return(y)

}
