\name{accessors}
\alias{accessors}
\alias{conditionSizes}
\alias{dispersions}
\alias{INICalls}
\alias{INIThreshold}
\alias{INIValues}
\alias{inputData}
\alias{logFC}
\alias{means}
\alias{normalizedData}
\alias{params}
\alias{posteriorProbs}
\alias{pvals}
\alias{responsibilities}
\alias{sampleNames}
\alias{sizeFactors}
\alias{sizeParameters}
\alias{transcriptNames}
\title{Accessors for a "DEXUSResult".}
\arguments{
  \item{object}{An instance of "DEXUSResult".}
}
\value{
  The accessor functions return a the matrices or vectors
  contained in the corresponding slot of the "DEXUSResult".
}
\description{
  These generic functions return the slots of an RNA-Seq
  analysis performed by DEXUS. The results of DEXUS are
  stored as an instance of \code{\link{DEXUSResult-class}}.
}
\examples{
data(dexus)
result <- dexus(countsBottomly[1:20,1:10])
transcriptNames(result)
sampleNames(result)
inputData(result)
normalizedData(result)
sizeFactors(result)
INIValues(result)
INIThreshold(result)
INICalls(result)
pvals(result)
responsibilities(result)
posteriorProbs(result)
logFC(result)
conditionSizes(result)
sizeParameters(result)
means(result)
dispersions(result)
params(result)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

