\name{getSizeNB}
\alias{getSizeNB}
\title{Maximum-likelihood and maximum-a-posteriori estimators
for the negative binomial distribution.}
\usage{
  getSizeNB(x, maxCyc = 1000, eta = 0, rmax = Inf,
    method = "bisection")
}
\arguments{
  \item{x}{The input data. Must be a numeric vector.}

  \item{maxCyc}{The maximum number of cycles of the numeric
  procedure to find the estimator. (Default = 1000).}

  \item{eta}{The weight of the exponential prior. The
  higher eta, the lower the estimate for the size
  parameter. Setting \env{eta} = 0 means that the prior is
  not used and, therefore, the maximum-likelihood estimator
  is calculated. (Default = 0).}

  \item{rmax}{Upper bound on the size parameter. This
  corresponds to a truncated exponential prior. If not used
  there is a non-zero probability that the estimator for
  the size parameter is \eqn{\infty}. (Default = Inf).}

  \item{method}{The procedure used to solve the equation

  \deqn{\sum_{k=1} ^N \psi (x_i+r) - N\psi(r)+N \log
  \left(\frac{r}{r+ 1/N \sum_{i=1}^N x_i} \right) - \eta =0
  }

  for \eqn{r}.

  This can either be "bisection" or "regula falsi".
  (Default="bisection").}
}
\value{
  "numeric" An estimate of the size parameter of the
  negative binomial distribution. The overdispersion
  parameter is the inverse of the size parameter of a
  negative binomial distribution
}
\description{
  Estimates the size parameter of a a negative binomial
  distribution from given data.
}
\details{
  Depending on the parameters you can either obtain the
  \emph{Maximum-likelihood estimator} or the
  \emph{maximum-a-posteriori estimator} using an
  exponential prior.

  \tabular{lc}{ maximum-likelihood estimator \tab \env{eta}
  = 0 \cr maximum-a-posteriori estimator \tab \env{eta} > 0
  \cr }

  By setting the variable \env{rmax} to a positive value
  one can enforce an upper bound on the parameter.

  The inverse of the size parameter is the overdispersion
  parameter.
}
\examples{
x <- rnbinom(mu=50, size=5, n=10)
getSizeNB(x)
}
\author{
  Guenter Klambauer \email{klambauer@bioinf.jku.at} and
  Thomas Unterthiner \email{unterthiner@bioinf.jku.at}
}

