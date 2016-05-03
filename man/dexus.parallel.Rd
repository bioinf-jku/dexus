\name{dexus.parallel}
\alias{dexus.parallel}
\title{A parallel version of DEXUS.}
\usage{
  dexus.parallel(X, ncores = 2, normalization = "RLE",
    ignoreIfAllCountsSmaller = 1, resultObject = "S4", ...)
}
\arguments{
  \item{X}{Either a vector of counts or a raw data matrix,
  where columns are interpreted as samples and rows as
  genomic regions.}

  \item{ncores}{The number of cores (CPUs) that will be
  used by the parallelization.}

  \item{normalization}{Normalization method to be used.
  (Default="RLE")}

  \item{ignoreIfAllCountsSmaller}{A transcript is
  considered as not expressed if all counts are smaller
  than the given value. (Default=1)}

  \item{resultObject}{Type of the result object; can either
  be a list ("list") or an instance of "DEXUSResult"
  ("S4"). (Default="S4").}

  \item{...}{Other options to be passed to \code{dexus()}.}
}
\value{
  "list"
}
\description{
  Speeds up DEXUS by using multiple processors. Uses the
  \code{parallel} package to parallelize a DEXUS call.
}
\examples{
data(dexus)
result <- dexus.parallel(countsPickrell[1:10, ],ncores=1)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

