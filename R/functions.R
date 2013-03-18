# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>


# makes sure the largest cluster will be assigned mainClassIdx
rename.largest.cluster <- function(clusters, n, mainClassIdx) {
	T <- table(factor(clusters,levels=1:n))
	x <- as.integer(names(T)[order(T,decreasing=TRUE)])
	y <- (mainClassIdx+(0:(n-1)))
	y[which(y>n)] <- y[which(y>n)]-n
	rclusters <- y[match(clusters,x)]
	return(rclusters)
}

# fast and numerically stable nbinom function
mydnbinom <- function(x,mu,size,log=FALSE){
	r <- size
	res <- lgamma(x+r)-lgamma(x+1)-lgamma(r)+ x*(log(mu)-log(mu+r)) + r*(log(r)-log(mu+r))
	if (log) return(res)
	else return(exp(res))
}

initialize.clusters.kmeans <- function(x, n, mainClassIdx, rmax, kmeansRuns=5) {
	# add noise to avoid "more cluster centers than distinct data points".
	k <- kmeans(log((x+0.01)/median(x+0.01))+abs(rnorm(n=length(x),sd=0.001)), n, nstart=kmeansRuns)
	clusters <- rename.largest.cluster(k$cluster, n, mainClassIdx)
	rmu <- sapply(1:n, function(i) {
				y = x[clusters==i];
				rr <- c(ifelse(length(y)>1, getSizeNB(y,rmax=rmax), 1.0),mean(y))
				return(rr)
			})
	
	r <- rmu[1, ]
	means <- pmax(rmu[2, ],0.01)
	
	if (length(which(table(clusters)==max(table(clusters))))!=1){
		r[c(mainClassIdx,(1:n)[-mainClassIdx])] <- r[order(means)]
		means[c(mainClassIdx,(1:n)[-mainClassIdx])] <- means[order(means)]
	}
	return(list(clusters=clusters,means=means,r=r))
}



initialize.clusters.quantiles <- function(x, n, mainClassIdx, rmax, dummy=NULL) {
	x <- x+rnorm(mean=0,sd=0.0001,n=length(x))
	r <- vector("numeric",n)
	mainQ <- 0.8
	
	t <- vector("numeric",n+1)
	if (n %% 2 ==0){
		skewness <- mean((x - mean(x))^3 )
		smallIntervals <- (1-mainQ)/(n-1)
		if (skewness >=0){
			t <- quantile(x, cumsum(c(rep(smallIntervals,n/2-1),mainQ,
											rep(smallIntervals,n/2))) )
		} else {
			t <- quantile(x, cumsum(c(rep(smallIntervals,n/2),mainQ,
											rep(smallIntervals,n/2-1))) )
		}
	} else {
		smallIntervals <- (1-mainQ)/(n-1)
		t <- quantile(x, cumsum(c(rep(smallIntervals,(n-1)/2),mainQ,
										rep(smallIntervals,(n-1)/2))) )
	}
	
	clusters <- apply(sapply(x,function(x) (x <= t)),2,
			function(x) length(which(x)))
	clusters <- rename.largest.cluster(clusters=clusters,n=n,
			mainClassIdx=mainClassIdx)
	r <- sapply(1:n, function(i) {
				y = x[clusters==i];
				ifelse(length(y)>1, getSizeNB(y,rmax=rmax), 
						y+abs(rnorm(n=1,sd=0.001)))
			})
	means <- sapply(1:n, function(i) {
				y = x[clusters==i];
				return(max(mean(y),1))
			})	
	return(list(clusters=clusters,means=means,r=r))
}


testMulticlass <- function(X, phi, design, design0 = NULL, normalization="none") {
	library(statmod)	
	library(stats)
	if (is.null(design0))
		design0 <- matrix(1, ncol(X))
	
	if (ncol(design0) > ncol(design)) {
		warning("H0 hypothesis design should be smaller than H1 design")
	}
	
	norm <- normalizeData(X, normalization)
	X <- norm$X
	
	res <- lapply(1:nrow(X), function(idx) {
				x <- X[idx, ]
				disp <- phi[idx]
				fit0 <- statmod::glmnb.fit(design0, x, disp, offset=log(norm$sizeFactors))			
				fit1 <- statmod::glmnb.fit(design, x, disp, offset=log(norm$sizeFactors))
				dev <- stats::deviance(fit0) - deviance(fit1)
				df <- length(fit1$coefficients) - length(fit0$coefficients)
				pval <- 1 - pchisq( dev, df)
				return(list(dev=dev, df=df, pval=pval, fit=fit1))
			})
	res
}



# basically estimateSizeFactorsForMatrix from DESeq
normalize.rle <- function(counts) {
	gm <- exp(rowMeans(log(counts)))
	f <- apply(counts, 2, function(col) median((col/gm)[gm > 0]))
	return(f)	
}


# upper-quartile normalization
normalize.uq <- function(counts, p=0.75) {
	c <- counts[rowSums(counts) > 0, ]
	uq <- apply(c,2,function(x) quantile(x,p=p))
	f <- uq/mean(uq) # scaling
	return(f)
}


