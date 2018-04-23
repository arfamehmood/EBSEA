\name{filterCounts}
\alias{filterCounts}
\title{ Filter Count Data
}
\description{
The exons are filtered based on their expression levels so that each exon has a mean greater than 1 or user defined mean.}
\usage{
filterCounts(x, mean)
}
\arguments{
\item{x}{ A numeric dataframe of counts in the sample with gene and
exon number as the row names and samples as the column names}
\item{mean}{ The mean value to filter the genes below the threshold}
}

\value{
{A dataframe of filtered counts of exons}
}
\seealso{\code{\link{EBSEA}}}
\examples{
data(origCounts)
res <- filterCounts(origCounts)
}
\keyword{ manip }
