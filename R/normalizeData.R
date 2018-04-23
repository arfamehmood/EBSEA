normalizeData <- function(x, group = NULL){
  message('Normalizing Data')
  norm.x <- DGEList(counts = x, group = group)
  norm.x <- calcNormFactors(norm.x, method = 'TMM')
  return(norm.x)
}
