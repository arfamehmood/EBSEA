\name{filterGenes}
\alias{filterGenes}
\title{Filter Gene List}
\description{The differentially expressed genes are filtered based on
the FC and FDR provided by the user.The default thresholds are FC => 1.25
and fdr <= 0.01}
\usage{
filterGenes(x, fc = 1.25, fdr = 0.01)
}
\arguments{
\item{x}{ A dataframe containing the gene statistics returned by EBSEA.}
\item{fc}{ A fold change threshold for the genes to be filtered. Default: 1.25.}
\item{fdr}{A FDR threshold for the genes to be filtered. Default: 0.01.}
}
\value{
{A list containing upregulated and downregulated genes.}
}
\seealso{\code{\link{EBSEA}}}
\examples{
data(origCounts)
group <- c('Group1', 'Group1', 'Group1', 'Group2', 'Group2', 'Group2', 'Group2')
result <- EBSEA(origCounts, group)
filteredGenes <- filterGenes(result$GeneTable)
}

\keyword{ manip }
