\name{sort}
\alias{sort}
\title{Sorting a DEXUS result.}
\arguments{
  \item{object}{An instance of "DEXUSResult".}
}
\value{
  An instance of "DEXUSResult".
}
\description{
  This function sorts the result object by I/NI values or
  p-values such that the transcripts with the highest I/NI
  value or the lowest p-value are ranked first.
}
\examples{
data(dexus)
res <- dexus(countsBottomly[1:100, ])
sort(res)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

