\name{INI}
\alias{INI}
\title{I/NI filtering of a DEXUS result.}
\arguments{
  \item{object}{An instance of "DEXUSResult".}

  \item{threshold}{A numeric determining the threshold for
  the I/NI values.}
}
\value{
  An instance of "DEXUSResult".
}
\description{
  This function filters the result object for informative
  transcripts. Transcripts with an I/NI value below the
  given threshold are filtered out.
}
\examples{
data(dexus)
res <- dexus(countsBottomly[1:100, ])
INI(res)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

