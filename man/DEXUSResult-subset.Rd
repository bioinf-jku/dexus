\name{`[`}
\alias{`[`}
\alias{`[`,DEXUSResult,character-method}
\alias{`[`,DEXUSResult,logical-method}
\alias{`[`,DEXUSResult,numeric-method}
\title{Subsetting a "DEXUSResult".}
\arguments{
  \item{x}{"DEXUSResult"}

  \item{i}{Either a numeric vector of indices or a
  character vector containing the transcript names.}
}
\value{
  An instance of "DEXUSResult".
}
\description{
  Information about specific transcripts can be accessed in
  the "DEXUSResult" object by using the standard brackets
  "[idx]" for subsetting. Either transcript names or
  transcript indices can be used.
}
\examples{
data(dexus)
res <- dexus(countsBottomly[1:100, ])
res["ENSMUSG00000000486"]
res[50:55]
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

