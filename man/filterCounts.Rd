\name{filterCounts}
\alias{filterCounts}
\title{ Filter Count Data
}
\description{
The exons are filtered based on their expression levels so that each exon has a
cpm (count per million) of more than 1 in user defined percent of the samples.}
\usage{
filterCounts(x, noOfSamples)
}
\arguments{
\item{x}{ A numeric dataframe of counts in the sample with gene and
exon number as the row names and samples as the column names}
\item{noOfSamples}{ Percentage of the number of samples
that should have cpm greater than 1.}
}

\value{
{A dataframe of filtered counts of exons}
}
\seealso{\code{\link{EBSEA}}}
\examples{
data(origCounts)
res <- filterCounts(origCounts, 20)
}
\keyword{ manip }
