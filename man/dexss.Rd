\name{dexss}
\alias{dexss}
\alias{DEXSS,}
\title{Detection of Differential Expression in a semi-supervised Setting}
\usage{
  dexss(X, nclasses = 2, G = 1, alphaInit, cyc = 20,
    labels, normalization = "RLE", kmeansIter = 10,
    ignoreIfAllCountsSmaller = 1, theta = 2.5, minMu = 0.5,
    rmax = 13, initialization = "kmeans",
    multiclassPhiPoolingFunction = NULL, quiet = FALSE,
    resultObject = "S4")
}
\arguments{
  \item{X}{either a vector of counts or a raw data matrix,
  where columns are interpreted as samples and rows as
  genomic regions. An instance of "countDataSet" is also
  accepted.}

  \item{nclasses}{The number of conditions, i.e. mixture
  components. (Default = 2)}

  \item{G}{The weight of the prior distribution of the
  mixture weights. Not used in the supervised case.
  (Default = 1).}

  \item{cyc}{Positive integer that sets the number of
  cycles of the EM algorithm. (Default = 20).}

  \item{alphaInit}{The initial estimates of the condition
  sizes, i.e., mixture weights. Not used in the supervised
  case. (Default = c(0.5,0.5)) .}

  \item{labels}{The labels for the classes, will be coerced
  into an integer. For this semi-supervised version the
  known labels/conditions must be coded as integers
  starting with 1. The samples with the label 1 will be
  considered as being in the "major condition". For the
  samples with unknown labels/conditions an "NA" must be
  set.}

  \item{normalization}{method used for normalizing the
  reads. "RLE" is the method used by (Anders and Huber,
  2010), "upperquartile" is the Upper-Quartile method by
  (Bullard et al., 2010), and none deactivates
  normalization. (Default = "RLE").}

  \item{kmeansIter}{number of times the K-Means algorithm
  is run. (Default = 10).}

  \item{ignoreIfAllCountsSmaller}{Ignores transcript for
  which all read counts are smaller than this value. These
  transcripts are considered as "not expressed" (Default =
  1).}

  \item{theta}{The weight of the prior on the size
  parameter or inverse dispersion parameter. Theta is
  adjusted to each transcript by dividing by the mean read
  count of the transcript. The higher \code{theta}, the
  lower \code{r} and the higher the overdispersion will be.
  (Default = 2.5).}

  \item{minMu}{Minimal mean for all negative binomial
  distributions. (Default = 0.5).}

  \item{rmax}{Maximal value for the size parameter. The
  inverse of this parameter is the lower bound on the
  dispersion. In analogy to (Anders and Huber, 2010) we use
  13 as default. (Default = 13).}

  \item{initialization}{Method used to find the initial
  clusters. Dexus can either use the quantiles of the
  readcounts of each gene or run k-means on the counts.
  (Default = "kmeans").}

  \item{multiclassPhiPoolingFunction}{In "multiClass" mode
  the dispersion is either estimated across all classes at
  once (NULL), or separately for each condition, i.e.,
  class. The size parameters or dispersion per class are
  then joined to one estimate by the mean ("mean"), minimum
  ("min") or maximum ("max"). In our investigations
  estimation across all classes at once performed best.
  (Default = NULL).}

  \item{quiet}{Logical that indicates whether dexus should
  report the steps of the algorithm. Supresses messages
  from the program if set to TRUE. (Default = FALSE).}

  \item{resultObject}{Type of the result object; can either
  be a list ("list") or an instance of "DEXUSResult"
  ("S4"). (Default="S4").}
}
\value{
  "list" or "DEXUSResult". A list containing the results
  and the parameters of the algorithm or an instance of
  "DEXUSResult".
}
\description{
  Performs the DEXSS algorithm for detection of
  differentially expressed genes in RNA-seq data for a
  semi-supervised setting, i.e. that the condition of some
  samples is known, and for some samples the condition is
  unkown.
}
\details{
  The read count \eqn{x} is explained by a finite mixture
  of negative binomials:

  \deqn{ p(x) = \sum_{i=1} ^n \alpha_i\ \mathrm{NB}(x;\
  \mu_i, r_i), }

  where \eqn{\alpha_i} is the weight of the mixture
  component, \eqn{\mathrm{NB}} is the negative binomial
  with mean parameter \eqn{\mu_i} and size parameter
  \eqn{r_i}. The parameters are selected by an EM algorithm
  in a Baysian framework.

  Each component in the mixture model corresponds to one
  condition.

  \itemize{ \item If the groups, conditions, replicate
  status, or labels are unknown, DEXUS tries to estimate
  these conditions. For each transcript DEXUS tries to
  explain the read counts by one negative binomial
  distribution. If this is possible, the transcript is
  explained by one condition and therefore it is not
  differentially expressed. If more than one negative
  binomial distribution is needed to explain the read
  counts of a transcript, this transcript indicates that it
  is differentially expressed. Evidence for differential
  expression is strong if a large amount of samples
  participate in each condition and the mean expression
  values are well separated. Both of these criteria are
  measured by the informative/non-informative (I/NI) call.

  \item If there are more than two groups given by the
  vector \env{labels}, DEXUS uses a generalized linear
  model to explain the data in analogy to (McCarthy, 2012).

  \item If there are two groups given by the vector
  \env{labels}, DEXUS uses the exact test for count data to
  test between the sample groups, as implemented by (Anders
  and Huber, 2010) in the package "DESeq". }
}
\examples{
data(dexus)
labels1 <- substr(colnames(countsBottomly),1,2)
labels2 <- c()
labels2[which(labels1=="D2")] <- 1
labels2[which(labels1=="B6")] <- 2
labels2[c(3,7,8,10,12,15)] <- NA
res <- dexss(countsBottomly[1:100, ],labels=labels2,nclasses=2,G=0)
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

  McCarthy, D. J., Chen, Y., and Smyth, G. K. (2012).
  \emph{Differential expression analysis of multifactor
  RNA-Seq experiments with respect to biological
  variation.} Nucleic Acids Res, 40(10), 4288-4297.
}

