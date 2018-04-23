\name{EBSEA}
\alias{EBSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Exon Based Startegy for Expression Analysis of genes
}
\description{
EBSEA takes as input unnormalized counts of exons, normalizes them and then
performs a two group comparison of the samples to detect differentially
expressed between the groups. Both paired or unpaired comparison are supported.
It calculates fold changes, p-values and false discovery rate of the genes
between the groups.
}
\usage{
EBSEA(data, group, paired = FALSE, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{data}{A dataframe of exon count data}
\item{group}{A vector indicating the sample groups in the experiment}
\item{paired}{A logical indicating whether the samples are paired or unpaired.
Default: FALSE}
\item{plot}{A logical indicating whether a volcano plot is visualized.
Default: FALSE}
}
\value{
{EBSEA returns a list of two  dataframes. ExonTable is a dataframe that
contains exon statistics including log fold change, p-values, adjusted p-values,
average expression and fold change. GeneTable is a dataframe that contains the
corresponding fold change, log fold change, p-values and false discovery rate.}
}


\references{
Laiho, A., & Elo, L. L. (2014). A note on an exon-based strategy to identify
differentially expressed genes in RNA-seq experiments. PloS One, 9(12), e115964.
}
\seealso{\code{\link{visualizeGenes}}}
\examples{
data(origCounts)
group <- c('Group1', 'Group1', 'Group1', 'Group2', 'Group2', 'Group2', 'Group2')
result <- EBSEA(origCounts, group)
}

\keyword{ distribution }
