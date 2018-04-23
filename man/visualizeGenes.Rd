\name{visualizeGenes}
\alias{visualizeGenes}
\title{Visualize Gene}
\description{
Plots for each exon of the gene entered by the user, the mean of the counts and the fold changes.
}
\usage{
visualizeGenes(gene, ebsea.out)
}

\arguments{
\item{gene}{Gene Name. The gene name should be the from the genes in
count data.}
\item{ebsea.out}{ Result object returned by EBSEA}
}
\value{
{A plot of mean counts and fold changes of exons of a gene.}
}
\seealso{\code{\link{EBSEA}}}
\examples{
data(origCounts)
group <- c('Group1', 'Group1', 'Group1', 'Group2', 'Group2', 'Group2', 'Group2')
result <- EBSEA(origCounts, group)
visualizeGenes('FBgn0000017', result)
}
\keyword{ device }

