\name{plot}
\alias{plot}
\title{Visualization of a result of the DEXUS algorithm.}
\arguments{
  \item{x}{An instance of "CNVDetectionResult"}

  \item{idx}{The indices or the transcript names of the
  transcripts that should be visualized as heatmap.}

  \item{cexSamples}{Size of the column labels, i.e. the
  samples.}

  \item{cexGenes}{Size of the row labels, i.e. the
  transcripts.}

  \item{newColNames}{renames the samples.}

  \item{type}{Mark the samples, that do not belong to the
  major class by crosses ("crosses"), or boxes ("boxes").}
}
\value{
  Generates a heatmap of the expression values of the
  top-ranked transcripts.
}
\description{
  Plots a heatmap of the log read counts of the top ranked
  genes.
}
\examples{
data(dexus)
r <- dexus(countsBottomly[1:100, ])
plot(r)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

