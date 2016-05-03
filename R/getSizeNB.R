# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>

#' @title Maximum-likelihood and maximum-a-posteriori estimators
#' for the negative binomial distribution.
#' 
#' @description Estimates the size parameter of a 
#' a negative binomial distribution from given data. 
#' 
#' @details Depending on the parameters you can
#' either obtain the \emph{Maximum-likelihood estimator} or the 
#' \emph{maximum-a-posteriori estimator} using an exponential prior.
#' 
#' \tabular{lc}{
#' maximum-likelihood estimator \tab \env{eta} = 0 \cr
#' maximum-a-posteriori estimator \tab \env{eta} > 0 \cr
#' }
#' 
#' By setting the variable \env{rmax} to a positive value one can enforce 
#' an upper bound on the parameter. 
#' 
#' The inverse of the size parameter is the overdispersion parameter.
#' 
#' @param x The input data. Must be a numeric vector. 
#' @param maxCyc The maximum number of cycles of the numeric procedure
#' to find the estimator. (Default = 1000).
#' @param eta The weight of the exponential prior. The higher eta, the lower
#' the estimate for the size parameter. Setting \env{eta} = 0 means
#' that the prior is not used and, therefore, the maximum-likelihood estimator
#' is calculated. 
#' (Default = 0).
#' @param rmax Upper bound on the size parameter. This corresponds to a truncated
#' exponential prior. If not used there is a non-zero probability that the
#' estimator for the size parameter is \eqn{\infty}. (Default = Inf).
#' @param method The procedure used to solve the equation
#' 
#' \deqn{\sum_{k=1} ^N \psi (x_i+r) - N\psi(r)+N \log \left(\frac{r}{r+
#'  1/N \sum_{i=1}^N x_i} \right) - \eta =0 }
#' 
#' for \eqn{r}.
#' 
#' This can either be "bisection" or "regula falsi". (Default="bisection").
#' 
#' 
#' @return "numeric" An estimate of the size parameter of the negative
#' binomial distribution. The overdispersion parameter is the inverse of
#' the size parameter of a negative binomial distribution
#' 
#' @examples 
#' x <- rnbinom(mu=50, size=5, n=10)
#' getSizeNB(x)
#' 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas 
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export

getSizeNB <- function(x, maxCyc=1000, eta=0, rmax=Inf, method="bisection"){
	if (length(x)==1){
		stop("Only one data point given. Cannot estimate variance!")
	}
	if(var(x)==0){
		return(min(Inf,rmax))
	} else if (var(x)/mean(x) <= 1) {
		return(min(Inf,rmax))
	} else {
		if (rmax==Inf) rmax <- -1
		if (method=="bisection") {
			res <- .C("find_r_bisection_ext",
					retval=as.numeric(0), 
					as.numeric(x),
					as.numeric(rep(1,length(x))),
					as.numeric(1),
					as.integer(length(x)),
					as.integer(1),
					as.integer(0),
					as.numeric(mean(x)),
					as.integer(0),
					as.integer(maxCyc),
					as.numeric(eta),
					as.numeric(rmax))$retval
		} else {
			if (rmax==Inf) rmax <- -1
			res <- .C("find_r_regulafalsi_ext",
					retval=as.numeric(0), 
					as.numeric(x),
					as.numeric(rep(1,length(x))),
					as.numeric(1),
					as.integer(length(x)),
					as.integer(1),
					as.integer(0),
					as.numeric(mean(x)),
					as.integer(0),
					as.integer(maxCyc),
					as.numeric(eta),
					as.numeric(rmax))$retval
		}	
		return(res)
	}
}

