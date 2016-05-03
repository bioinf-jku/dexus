\name{INIThreshold<-}
\alias{INIThreshold<-}
\alias{INIThreshold-set}
\title{Set the I/NI threshold.}
\arguments{
  \item{object}{An instance of "DEXUSResult".}

  \item{value}{A numeric to be used for thresholding the
  I/NI values.}
}
\value{
  \code{INIThreshold<-} returns an instance of
  "DEXUSResult".
}
\description{
  This generic function sets the threshold of the I/NI
  value. Transcripts with I/NI values above the I/NI
  threshold are considered as differentially expressed. The
  results of DEXUS are stored as an instance of
  \code{\link{DEXUSResult-class}}.
}
\examples{
data(dexus)
result <- dexus(countsBottomly[1:20,1:10])
INIThreshold(result) <- 0.1
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

