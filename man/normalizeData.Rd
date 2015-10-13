\name{normalizeData}
\alias{normalizeData}
\title{Normalization of RNA-Seq count data.}
\usage{
  normalizeData(X, normalization)
}
\arguments{
  \item{X}{data a raw data matrix, where' columns are
  interpreted as samples and rows as genomic regions.}

  \item{normalization}{method used for normalizing the
  reads. RLE is the method used by (Anders and Huber,
  2010), upperquartile is the Upper-Quartile method from
  (Bullard et al., 2010), and none deactivates
  normalization. (Default = "RLE").}
}
\value{
  "list" A list containing the normalized data (in its "X"
  component) as well as the size-factors used for the
  normalization ("sizeFactors").
}
\description{
  Normalizes RNA-seq count data using previously published
  approaches. Each samples' read counts are corrected by a
  normalizing factor. The options are "RLE" by (Anders and
  Huber, 2010), and "upperquartile" by (Bullard et al.,
  2010).
}
\examples{
data(dexus)
norm <- normalizeData(countsBottomly,"RLE")
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}
\references{
  Anders, S. and Huber, W. (2010). \emph{Differential
  expression analysis for sequence count data.} Genome
  Biol, 11(10), R106.

  Bullard, J. H., Purdom, E., Hansen, K. D., and Dudoit, S.
  (2010). \emph{Evaluation of statistical methods for
  normalization and differential expression in mRNA-seq
  experiments.} BMC Bioinformatics, 11, 94.
}

