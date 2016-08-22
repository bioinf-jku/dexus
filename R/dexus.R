# Copyright (C) 2013 Klambauer Guenter and Thomas Unterthiner
# <klambauer@bioinf.jku.at>

#' @title A parallel version of DEXUS.
#'
#' @description Speeds up DEXUS by using multiple processors.
#' Uses the \code{parallel} package to parallelize a DEXUS call.
#' 
#' @param X Either a vector of counts or a raw data matrix, where
#'	columns are interpreted as samples and rows as genomic regions.
#' @param ncores The number of cores (CPUs) that will be used by 
#' the parallelization.
#' @param normalization Normalization method to be used. (Default="RLE")
#' @param ignoreIfAllCountsSmaller A transcript is considered as not expressed
#' if all counts are smaller than the given value. (Default=1)
#' @param resultObject Type of the result object; can either be a list ("list")
#' or an instance of "DEXUSResult" ("S4"). (Default="S4").
#' @param ... Other options to be passed to \code{dexus()}.
#' @return "list" 
#' 
#' @examples 
#' data(dexus)
#' result <- dexus.parallel(countsPickrell[1:10, ],ncores=1)
#' 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas 
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export


dexus.parallel <- function(X, ncores=2,
		normalization = "RLE",
		ignoreIfAllCountsSmaller=1,resultObject="S4", ...){
	library(parallel)
	### check input ############################################################
	if (is.matrix(X)){
		if (is.null(rownames(X))){
			rownames(X) <- paste("Transcript_",1:nrow(X),sep="")
		} else {
			if (length(unique(rownames(X)))!=length(rownames(X)))
				stop("Transcript names must be unique.")
		}
		if (is.null(colnames(X))){
			colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		} else {
			if (length(unique(colnames(X)))!=length(colnames(X)))
				stop("Sample names must be unique.")
		}
	} else if (is.vector(X)) {
		X <- matrix(X,nrow=1)
		rownames(X) <- "Transcript1"
		colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		normalization <- "none"
	} else if (inherits(X,"CountDataSet")) {
		library(DESeq)
		object <- X
		X <- DESeq::counts(X)
		labels <- DESeq::conditions(object)
		if (is.null(rownames(X))){
			rownames(X) <- paste("Transcript_",1:nrow(X),sep="")
		} else {
			if (length(unique(rownames(X)))!=length(rownames(X)))
				stop("Transcript names must be unique.")
		}
		if (is.null(colnames(X))){
			colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		} else {
			if (length(unique(colnames(X)))!=length(colnames(X)))
				stop("Sample names must be unique.")
		}
	} else {
		stop("Input data types must be a matrix or vector!")
	}
	
	if (!(normalization %in% c("none","RLE","upperquartile"))){
		stop("\"normalization\" must be \"none\",\"RLE\" or \"upperquartile\".")
	}
	###
	
	if (ncores==1 | !require(parallel, quietly=TRUE)) {
		message("Switching to single-threaded code.")
		return(dexus(X, normalization=normalization, 
						ignoreIfAllCountsSmaller=ignoreIfAllCountsSmaller,
						resultObject=resultObject,...))
	}
	
	trows <- nrow(X)
	N <- ncol(X)
	X.raw <- X
	norm.method  <-  normalization
	
	norm <- normalizeData(X, norm.method)
	X <- norm$X
	
	exclIdx <- (apply(X,1,max) < ignoreIfAllCountsSmaller)
	if (length(exclIdx)) 
		message(paste("Filtered out ", 100*round(mean(exclIdx),4),
						"% of the genes due to low counts"))
	
	d <- as.integer(nrow(X) / ncores)
	res.par <- mclapply(1:ncores, function(i) {			
				a <- (i-1)*d+1
				b <- i*d
				if (i == ncores)
					b <- nrow(X)
				return (dexus(X[a:b, ], normalization="none",
								ignoreIfAllCountsSmaller=ignoreIfAllCountsSmaller, 
								quiet=TRUE,resultObject="list", ...))
			}, mc.cores=ncores)	
	
	res <- list(p=NULL, alpha=NULL, A=NULL, resp=NULL, r=NULL, means=NULL,
			sds=NULL,
			logfc=NULL, INI=NULL, pval=NULL, X=X.raw, X.norm=X, 
			sizeFactors=norm$sizeFactors, rINIT=NULL, meansINIT=NULL)
	for (i in 1:ncores){
		res$p = rbind(res$p, res.par[[i]]$p)
		res$alpha = rbind(res$alpha, res.par[[i]]$alpha)
		res$resp = rbind(res$resp, res.par[[i]]$resp)
		res$r = rbind(res$r, res.par[[i]]$r)
		res$means = rbind(res$means, res.par[[i]]$means)
		res$sds = rbind(res$sds, res.par[[i]]$sds)
		res$logfc = rbind(res$logfc, res.par[[i]]$logfc)
		res$INI = c(res$INI, res.par[[i]]$INI)
		res$pval = c(res$pval, res.par[[i]]$pval)
		res$rINIT = rbind(res$rINIT, res.par[[i]]$rINIT)
		res$meansINIT = rbind(res$meansINIT, res.par[[i]]$meansINIT)
	}
	if (!is.null(res.par[[1]]$A))
		res$A <- array(unlist(sapply(1:ncores, function(i) res.par[[i]]$A)), 
				dim=c(ncol(res$alpha),N,trows))
	else
		res$A <- NULL
	res$params <- res.par[[1]]$params
	res$params$normalization <- normalization
	res$params$resultObject <- resultObject
	
	if (resultObject!="S4"){
		res <- res[ c("params","X","X.norm", "sizeFactors",
						"A","resp","logfc","INI","pval",
						"rINIT", "meansINIT","p","alpha",
						"r","means","sds")]
		return(res)
	} else {
		res2 <- new("DEXUSResult")
		res2@transcriptNames	 <- rownames(X)
		res2@sampleNames		 <- colnames(X)
		res2@inputData		 <- res$X
		res2@normalizedData   <- res$X.norm  
		res2@sizeFactors		 <- norm$sizeFactors
		res2@INIValues        <- res$INI
		res2@INICalls		 <- res$INI>0.1
		res2@INIThreshold 	 <- 0.1
		res2@pvals			 <- res$pval
		res2@responsibilities <- res$resp
		res2@posteriorProbs	 <- res$A	
		res2@logFC			 <- res$logfc 	
		res2@conditionSizes	 <- res$alpha	
		res2@sizeParameters   <- res$r		
		res2@means            <- res$means			
		res2@dispersions	     <- (1/res$r)		
		res2@params           <- res$params 
		return(res2)
	}
}



#' @title Detection of Differential Expression in an Unsupervised Setting

#' @description Performs the DEXUS algorithm for detection of differentially 
#' expressed genes in RNA-seq data for a) unknown conditions, b) multiple 
#' known conditions, and c) two known conditions.
#' 
#' @details The read count \eqn{x} is explained by
#' a finite mixture of negative binomials:
#' 
#' \deqn{ 
#' p(x) = \sum_{i=1} ^n \alpha_i\  \mathrm{NB}(x;\ \mu_i, r_i),
#' }
#' 
#' where \eqn{\alpha_i} is the weight of the mixture component, \eqn{\mathrm{NB}}
#' is the negative binomial with mean parameter \eqn{\mu_i} and size parameter
#' \eqn{r_i}. The parameters are selected by an EM algorithm in a Baysian 
#' framework.
#' 
#' Each component in the mixture model corresponds to one condition.
#' 
#' \itemize{
#' \item If the groups, conditions, replicate status, or labels are unknown, DEXUS
#' tries to estimate these conditions. For each transcript DEXUS tries to
#' explain the read counts by one negative binomial distribution. If this is
#' possible, the transcript is explained by one condition and therefore it 
#' is not differentially expressed. If more than one negative binomial 
#' distribution is needed to explain the read counts of a transcript, this
#' transcript indicates that it is differentially expressed. Evidence for
#' differential expression is strong if a large amount of samples participate
#' in each condition and the mean expression values are well separated. Both
#' of these criteria are measured by the informative/non-informative (I/NI)
#' call. 
#' 
#' \item If there are more than two groups given by the vector \env{labels}, 
#' DEXUS uses a generalized linear model to explain the data in analogy to
#' (McCarthy, 2012). 
#' 
#' \item If there are two groups given by the vector \env{labels}, DEXUS
#' uses the exact test for count data to test between the 
#' sample groups, as implemented by (Anders and Huber, 2010) in the package 
#' "DESeq".
#' }
#' 
#'  
#' 
#' @param X either a vector of counts or a raw data matrix, where
#' columns are interpreted as samples and rows as genomic regions. An instance
#' of "countDataSet" is also accepted.
#' @param nclasses The number of conditions, i.e. mixture components. (Default = 2)
#' @param alphaInit The initial estimates of the condition sizes, i.e., mixture weights.
#' Not used in the  supervised case. (Default = c(0.5,0.5)) .
#' @param G The weight of the prior distribution of the mixture weights. 
#' Not used in the supervised case. (Default = 1).
#' @param cyc Positive integer that sets the number of cycles of the EM algorithm. 
#' (Default =  20).
#' @param labels labels for the classes, will be coerced into a factor
#' by \code{as.factor}. Can either be a factor, character or integer. If
#' this vector is given, supervised detection is used. If this vector is set 
#' to NULL the unsupervised detection is performed. (Default=NULL).
#' @param normalization method used for normalizing the reads. 
#' "RLE" is the method used by (Anders and Huber, 2010),
#' "upperquartile" is the Upper-Quartile method by (Bullard et al., 2010), and none
#' deactivates normalization. (Default = "RLE").
#' @param kmeansIter number of times the K-Means algorithm is run. (Default = 10).
#' @param ignoreIfAllCountsSmaller Ignores transcript for which all read counts
#' are smaller than this value. These transcripts are considered as "not expressed" 
#' (Default = 1).
#' @param theta The weight of the prior on the size parameter or inverse
#' dispersion parameter. Theta is adjusted to each transcript by dividing
#' by the mean read count of the transcript. The higher \code{theta}, the lower \code{r} and
#' the higher the overdispersion will be. (Default = 2.5).
#' @param minMu Minimal mean for all negative binomial distributions. (Default = 0.5).
#' @param rmax Maximal value for the size parameter. The inverse of this parameter
#' is the lower bound on the dispersion. In analogy to (Anders and Huber, 2010) we
#' use 13 as default. (Default = 13).
#' @param initialization Method used to find the initial clusters. 
#' Dexus can either use the quantiles of the readcounts of each gene or 
#' run k-means on the counts. (Default = "kmeans").
#' @param multiclassPhiPoolingFunction In "multiClass" mode the dispersion 
#' is either estimated across all classes at once (NULL), or separately for 
#' each condition, i.e., class. The size parameters or dispersion per class are
#' then joined to one estimate by the mean ("mean"), minimum ("min") or maximum ("max").
#' In our investigations estimation across all classes at once performed best.
#' (Default = NULL). 
#' @param quiet Logical that indicates whether dexus should report the steps
#' of the algorithm. Supresses messages from the program if set to TRUE. (Default = FALSE). 
#' @param resultObject Type of the result object; can either be a list ("list")
#' or an instance of "DEXUSResult" ("S4"). (Default="S4").
#' @return "list" or "DEXUSResult".  
#' A list containing the results and the parameters of the algorithm or an
#' instance of "DEXUSResult".
#' @references 
#' 
#' Anders, S. and Huber, W. (2010). \emph{Differential expression analysis
#' for sequence count data.} Genome Biol, 11(10), R106.
#' 
#' Bullard, J. H., Purdom, E., Hansen, K. D., and Dudoit, S. (2010).
#' \emph{Evaluation of statistical methods for normalization and differential
#'  expression in mRNA-seq experiments.} 
#' BMC Bioinformatics, 11, 94.
#' 
#' McCarthy, D. J., Chen, Y., and Smyth, G. K. (2012). 
#' \emph{Differential expression analysis of multifactor RNA-Seq experiments
#' with respect to biological variation.} Nucleic Acids Res, 40(10),
#' 4288-4297.
#' 
#' @examples 
#' data(dexus)
#' result <- dexus(countsMontgomery[1:10, ])
#' 
#' @aliases DEXUS, dexus
#' @useDynLib dexus
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas 
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export

dexus <- function(X, nclasses=2, alphaInit, G=1, 
		cyc=20, labels=NULL, normalization = "RLE", 
		kmeansIter=10, ignoreIfAllCountsSmaller=1, theta=2.5, minMu=0.5,rmax=13.0,
		initialization="kmeans", 
		multiclassPhiPoolingFunction=NULL, quiet=FALSE, resultObject="S4"){
	
	### check input ############################################################
	if (is.matrix(X)){
		if (is.null(rownames(X))){
			rownames(X) <- paste("Transcript_",1:nrow(X),sep="")
		} else {
			if (length(unique(rownames(X)))!=length(rownames(X)))
				stop("Transcript names must be unique.")
		}
		if (is.null(colnames(X))){
			colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		} else {
			if (length(unique(colnames(X)))!=length(colnames(X)))
				stop("Sample names must be unique.")
		}
	} else if (is.vector(X)) {
		X <- matrix(X,nrow=1)
		rownames(X) <- "Transcript1"
		colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		normalization <- "none"
	} else if (inherits(X,"CountDataSet")) {
		library(DESeq)
		object <- X
		X <- DESeq::counts(X)
		labels <- DESeq::conditions(object)
		if (is.null(rownames(X))){
			rownames(X) <- paste("Transcript_",1:nrow(X),sep="")
		} else {
			if (length(unique(rownames(X)))!=length(rownames(X)))
				stop("Transcript names must be unique.")
		}
		if (is.null(colnames(X))){
			colnames(X) <- paste("Sample_",1:ncol(X),sep="")
		} else {
			if (length(unique(colnames(X)))!=length(colnames(X)))
				stop("Sample names must be unique.")
		}
	} else {
		stop("Input data types must be a matrix or vector!")
	}
	
	
	if (!is.null(labels)){
		if (is.integer(labels) | is.numeric(labels)){
			labels <- as.integer(labels-1)
		} else if (is.factor(labels)){
			labels <- as.integer(labels)-1
		} else if (is.character(labels)){
			labels <- as.integer(as.factor(labels))-1
		} else {
			stop("Labels must be integer, factor or character.\n")
		}
		if (length(unique(labels))!=nclasses){
			message("Resetting number of classes due to labels.")
			nclasses <- length(unique(labels))
		}
	}
	

	
	if (!is.numeric(nclasses) | nclasses < 1 | nclasses!=as.integer(nclasses)){
		stop("\"nclasses\" must be integer greater or equal to 1.")
	}
	if (nclasses > ncol(X))
		stop("Number of conditions greater than the number of observations.")
	
	if (missing(alphaInit))
		alphaInit <- rep(1,nclasses)
	if (!is.numeric(alphaInit) | any(alphaInit < 0) | length(alphaInit)!=nclasses){
		message("\"alphaINIT\" must be numeric greater and all values greater 0.")
		stop("Length of this vector must be \"nclasses\"")
	}
	
	if (!is.numeric(G) | length(G) !=1 | G < 0){
		stop("\"G\" must be a numeric greater or equal to 0.")
	}
	if (!is.numeric(cyc) | length(cyc) !=1 | cyc < 1){
		stop("\"cyc\" must be an integer greater or equal to 1.")
	}
	
	
	
	if (!(normalization %in% c("none","RLE","upperquartile"))){
		stop("\"normalization\" must be \"none\",\"RLE\" or \"upperquartile\".")
	}
	if (!is.numeric(kmeansIter) | length(kmeansIter) !=1 | kmeansIter < 1){
		stop("\"kmeansIter\" must be an integer greater or equal to 1.")
	}
	if (!is.numeric(ignoreIfAllCountsSmaller) | length(ignoreIfAllCountsSmaller) !=1){
		stop("\"ignoreIfAllCountsSmaller\" must be a numeric value.")
	}
	if (!is.numeric(theta) | length(theta) !=1 | theta < 0){
		stop("\"theta\" must be a positive numeric value.")
	}
	if (!is.numeric(minMu) | length(minMu) !=1 | minMu < 0){
		stop("\"minMu\" must be a positive numeric value.")
	}
	if (!is.numeric(rmax) | length(rmax) !=1){
		stop("\"rmax\" must be a positive numeric value.")
	}
	if (rmax==Inf) rmax <- -1
	if (!(initialization %in% c("kmeans","quantiles"))){
		stop("\"initialization\" must be \"kmeans\" or \"quantiles\".")
	}
	if (!(is.null(multiclassPhiPoolingFunction))){
		if (!(multiclassPhiPoolingFunction %in% c("min","max","mean"))){
			stop(paste("\"multiclassPhiPoolingFunction\" must be \"min\"",
							"\"max\", or \"mean\" or NULL"))
		}
	}
	if (!is.logical(quiet)){
		stop("\"quiet\" must be a logical.")
	}
	
	
	
	############################################################################
	
	
	# Normalization
	X.raw <- X
	norm <- normalizeData(X, normalization)
	X <- norm$X
	
	gamma <- rep(1,nclasses)
	gamma[1] <- 1+G
	#eps1 <- 10^(-(1:nclasses)-6)
	#gamma <- gamma + eps1
	mainClassIdx <- 1
	trows <- nrow(X)
	exclIdx <- (apply(X,1,max) < ignoreIfAllCountsSmaller)
	
	if (all(exclIdx)){
		message("All transcripts filtered out because of low counts.")
		stop("Check your count matrix or change \"ignoreIfAllCountsSmaller\".")
	}
	
	if (length(exclIdx) & !quiet) 
		message(paste("Filtered out",100*round(mean(exclIdx),4),
						"% of the genes due to low counts"))
	X.all <- X
	X <- pmax(X,minMu)
	X <- X[!exclIdx, ,drop=FALSE]
	N <- ncol(X)
	g <- nrow(X)
	etaPerGene <- theta/((1+rowMeans(X)))
	varToMeanT <- 1
	
	if (is.null(labels)){
		### Unsupervised############################################################
		mode <- "unsupervised"
		
		if (!quiet) message("Unsupervised mode.")
		clusterFunc <- switch(initialization,
				kmeans=initialize.clusters.kmeans,
				quantiles=initialize.clusters.quantiles)
		
		
		clusterResults <- (apply(X,1,clusterFunc, n=nclasses, 
							mainClassIdx=mainClassIdx, rmax=rmax,kmeansIter))
		
		meansINIT <- sapply(clusterResults, function(x) x$means)
		rINIT <- sapply(clusterResults, function(x) x$r)
		alphaInit <- as.numeric(alphaInit/sum(alphaInit))
		
		
		res <- .Call("dexus",as.vector(t(X)), nrow(X), ncol(X),
				alphaInit, as.vector(rINIT), as.vector(meansINIT),
				as.numeric(gamma), as.integer(nclasses), as.integer(cyc),
				as.numeric(varToMeanT), as.numeric(etaPerGene),
				as.numeric(minMu),as.numeric(rmax))	
		
		## rerun if the major class is not mainClassIdx after learning
		alpha <- matrix(res$alpha,nrow=g,byrow=TRUE)
		A <- array(res$A,dim=c(nclasses,N,g))  ##right!
		r <- matrix(res$r,nrow=g,byrow=TRUE)
		means <-  matrix(res$means,nrow=g,byrow=TRUE)
		resp <- t(sapply(1:g,function(i) apply(A[,,i],2,which.max) ))	
		rerunIdx <- which(apply(resp,1,function(x) any(table(x) > table(x)[mainClassIdx])))
		if (length(rerunIdx)>0){
			reorder <- t(apply(resp,1,function(x){
								h <- hist(x,breaks=seq(0.5,nclasses+0.5,1),plot=FALSE)$counts
								order(h,decreasing=TRUE)[order(gamma,decreasing=TRUE)]
							}
					))
			res2 <- .Call("dexus",as.vector(t(X[rerunIdx, ,drop=FALSE])), 
					length(rerunIdx), ncol(X),
					alphaInit, 
					as.vector(t(r[rerunIdx,reorder[rerunIdx, ,drop=FALSE]])), 
					as.vector(t(means[rerunIdx,reorder[rerunIdx, ,drop=FALSE] ])),
					as.numeric(gamma), as.integer(nclasses), as.integer(cyc),
					as.numeric(varToMeanT), as.numeric(etaPerGene),
					as.numeric(minMu),as.numeric(rmax))
			
			#implant res2 into res... 
			implantIdx <- as.vector(sapply(rerunIdx*nclasses,
							function(x) (x-nclasses+1):x ))
			res$r[implantIdx] <- res2$r
			res$means[implantIdx] <- res2$means
			res$alpha[implantIdx] <- res2$alpha
			
			implantIdx2 <- as.vector(sapply((rerunIdx-1)*N*nclasses,
							function(x) (x+1):(x+N*(nclasses) )))
			res$A[implantIdx2] <- res2$A
			#check: res$A[implantIdx2]==as.vector(A[,,rerunIdx])
		}	
		
		pval <- rep(NA,nrow(X))
	} else {
		###   Supervised########################################################
		
		if (!quiet) message("Supervised mode.")
		
		
		if (nclasses > 2) {  
			library(statmod)
			library(stats)
			######### supervised multiclass ####################################
			mode <- "multi-class"
			
			design <- stats::model.matrix(~as.factor(labels))
			
			if (is.null(multiclassPhiPoolingFunction)) {
				# one dispersion for all components
				#nclasses <- 1
				alphaInit <- rep(1,1)
				alphaInit <- as.numeric(alphaInit/sum(alphaInit))
				meansINIT <- t(matrix(1, nrow=g, ncol=1))
				rINIT <- t(matrix(1,nrow=g,ncol=1))
				
				res <- .Call("dexusS",as.vector(t(X)), nrow(X), ncol(X),
						alphaInit, as.numeric(rINIT), as.numeric(meansINIT),
						as.numeric(gamma),
						as.integer(1), as.integer(1),as.numeric(varToMeanT),
						as.integer(labels), as.numeric(etaPerGene),as.numeric(minMu),
						as.numeric(rmax))
				r <- res$r
				phi <- 1/r
				res <- .Call("dexusS",as.vector(t(X)), nrow(X), ncol(X),
						alphaInit, as.numeric(rINIT), as.numeric(meansINIT),
						as.numeric(gamma),
						as.integer(nclasses), as.integer(1),as.numeric(varToMeanT),
						as.integer(labels), as.numeric(etaPerGene),as.numeric(minMu),
						as.numeric(rmax))	
				res$r <- as.vector(sapply(r,rep,nclasses))
				
			} else {
				# different phi for each component
				alphaInit <- rep(1,nclasses)
				alphaInit <- as.numeric(alphaInit/sum(alphaInit))	
				meansINIT <- t(matrix(1, nrow=g, ncol=nclasses))
				rINIT <- t(matrix(1,nrow=g,ncol=nclasses))
				
				res <- .Call("dexusS",as.vector(t(X)), nrow(X), ncol(X),
						alphaInit, as.numeric(rINIT), as.numeric(meansINIT),
						as.numeric(gamma),
						as.integer(nclasses), as.integer(1),as.numeric(varToMeanT),
						as.integer(labels), as.numeric(etaPerGene),as.numeric(minMu),
						as.numeric(rmax))	
				phi <- 1/(matrix(res$r,nrow=g,byrow=TRUE))
				phi <- apply(phi, 1, multiclassPhiPoolingFunction)
			}
			
			tmp <- testMulticlass(X, phi, design, normalization="none")
			pval <- sapply(tmp, function(x) x$pval)
			
		}  else if (nclasses == 2) {
			######### two classes class ########################################
			mode <- "two-class"
			
			meansINIT <- t(matrix(1, nrow=g, ncol=nclasses))
			rINIT <- t(matrix(1,nrow=g,ncol=nclasses))
			alphaInit <- rep(1,nclasses)
			alphaInit <- as.numeric(alphaInit/sum(alphaInit))	
			res <- .Call("dexusS",as.vector(t(X)), nrow(X), ncol(X),
					alphaInit, as.numeric(rINIT), as.numeric(meansINIT),as.numeric(gamma),
					as.integer(nclasses), as.integer(1),as.numeric(varToMeanT),
					as.integer(labels), as.numeric(etaPerGene),as.numeric(minMu),
					as.numeric(rmax))	
			
			
			if (!quiet) message("Performing calculation of p values.")
			idxGroupA <- which(labels==0)
			idxGroupB <- which(labels==1)
			if (length(idxGroupA)==0 | length(idxGroupB)==0) browser()
			
			# if rmax was set to "Inf"
			poisIdx <- which(res$r < 0) 
			res$r[poisIdx] <- Inf
			r <- matrix(res$r,nrow=g,byrow=TRUE)
			
			dispsA <- 1/r[,1]+1e-15
			dispsB <- 1/r[,2]+1e-15
			
			pval <- nbinomTestForMatrices(X.raw[!exclIdx,idxGroupA,drop=FALSE], X.raw[!exclIdx,idxGroupB,drop=FALSE],
					norm$sizeFactors[idxGroupA], norm$sizeFactors[idxGroupB], dispsA, dispsB )
						
			# Readjusting the low-ranked pvalues
			idxLR <- which(pval >= 0.99 | is.na(pval)) #LR:= Low Ranked
			if (length(idxLR)>0){
				pvalLR <- sort(runif(n=length(idxLR),0.01,1.00))
				rankO <- order(abs(rowSums(X[idxLR,idxGroupA,drop=FALSE])-
										rowSums(X[idxLR,idxGroupB,drop=FALSE])),decreasing=TRUE)
				idxLR <- idxLR[rankO]
				pval[idxLR] <- pvalLR 
			}
			
			
		} else {
			message("Labels must indicate at least two classes!")
			stop("For unsupervised mode set labels to NULL!")
		}
	}
	
	# extracting from result object
	alpha <- matrix(res$alpha,nrow=g,byrow=TRUE)
	if (mode=="unsupervised" & nclasses >2){
		oo <- t(apply(alpha,1,order,decreasing=TRUE))
		alpha <- t(sapply(seq(g), function(x) alpha[x, oo[x, ]] ))
	}
	
	A <- array(res$A,dim=c(nclasses,N,g))  ##right!
	if (mode=="unsupervised" & nclasses >2){
		A <- lapply(seq(g), function(x) A[oo[x, ], ,x] )
		A <- array(unlist(A),dim=c(nclasses,N,g))
	}
	
	
	poisIdx <- which(res$r < 0) 
	res$r[poisIdx] <- Inf
	r <- matrix(res$r,nrow=g,byrow=TRUE)
	if (mode=="unsupervised" & nclasses >2)
		r <- t(sapply(seq(g), function(x) r[x, oo[x, ]] ))
	
	means <-  matrix(res$means,nrow=g,byrow=TRUE)
	if (mode=="unsupervised" & nclasses >2)
		means <- t(sapply(seq(g), function(x) means[x, oo[x, ]] ))
	
	p <- res$means / (res$means + res$r)
	sds <- sqrt(res$means * (res$means + res$r) / res$r)
	sds[which(r==Inf)] <- sqrt(res$means[which(r==Inf)])
	sds <- matrix(sds,nrow=g,byrow=TRUE)
	resp <- t(sapply(1:g,function(i) apply(A[,,i],2,which.max) ))
	fc <- log((means+1e-4)/(means[,mainClassIdx]+1e-4))
	INI <- sapply(1:g, function(i) sum(alpha[i, ] * abs(fc)[i, ]))
	
	### reorder 
	
	
	### re-insert excluded genes
	pFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	pFinal[exclIdx, ] <- rep(NA,nclasses)
	pFinal[!exclIdx, ] <- p
	rownames(pFinal) <- rownames(X.all)
	colnames(pFinal) <- paste("Condition",1:nclasses,sep="_")
	
	
	alphaFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	alphaFinal[exclIdx, ] <- rep(0,nclasses)
	alphaFinal[exclIdx,mainClassIdx] <- 1
	alphaFinal[!exclIdx, ] <- alpha
	rownames(alphaFinal) <- rownames(X.all)
	colnames(alphaFinal) <- paste("Condition",1:nclasses,sep="_")
	
	
	AFinal <- array(NA,dim=c(nclasses,N,trows))
	AFinal[, ,exclIdx] <- matrix(NA,nrow=nclasses,ncol=N)
	AFinal[, ,!exclIdx] <- A
	
	respFinal <- matrix(NA,nrow=trows,ncol=N)
	respFinal[exclIdx, ] <- rep(mainClassIdx,N)
	respFinal[!exclIdx, ] <- resp
	rownames(respFinal) <- rownames(X.all)
	colnames(respFinal) <- colnames(X.all)
	
	
	rFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	rFinal[exclIdx, ] <- rep(NA,nclasses)
	rFinal[!exclIdx, ] <- r
	rownames(rFinal) <- rownames(X.all)
	colnames(rFinal) <- paste("Condition",1:nclasses,sep="_")
	
	meansFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	if (any(exclIdx))
		meansFinal[exclIdx, ] <- 
				t(sapply(rowMeans(X.all[exclIdx, ,drop=FALSE]),rep,nclasses))
	meansFinal[!exclIdx, ] <- means
	rownames(meansFinal) <- rownames(X.all)
	colnames(meansFinal) <- paste("Condition",1:nclasses,sep="_")
	
	sdsFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	sdsFinal[exclIdx, ] <- rep(NA,nclasses)
	sdsFinal[!exclIdx, ] <- sds
	rownames(sdsFinal) <- rownames(X.all)
	colnames(sdsFinal) <- paste("Condition",1:nclasses,sep="_")
	
	fcFinal <- matrix(NA,nrow=trows,ncol=nclasses)
	fcFinal[exclIdx, ] <- rep(NA,nclasses)
	fcFinal[!exclIdx, ] <- fc
	rownames(fcFinal) <- rownames(X.all)
	colnames(fcFinal) <- paste("Condition",1:nclasses,sep="_")
	
	INIFinal <- vector("numeric", length=trows)
	INIFinal[!exclIdx] <- INI
	if (length(which(exclIdx))>0)
		INIFinal[exclIdx] <- runif(n=length(which(exclIdx)),min=0,max=min(INI))
	names(INIFinal) <- rownames(X.all)
	
	pvalFinal <- rep(1, length=trows)
	pvalFinal[!exclIdx] <- pval
	if (all(is.na(pval))){
		pvalFinal[exclIdx] <- NA
	} else{
		pvalFinal[exclIdx] <- runif(n=length(which(exclIdx)),0.01,1.00)
	}
	names(pvalFinal) <- rownames(X.all)
	
	
	params <- list(nclasses, alphaInit, G, cyc, labels, normalization,
			kmeansIter, ignoreIfAllCountsSmaller, theta, minMu, rmax, initialization,
			mode, multiclassPhiPoolingFunction,quiet,resultObject)
	names(params) <- c("nclasses","alphaINIT","G","cyc","labels",
			"normalization","kmeansIter","ignoreIfAllCountsSmaller",
			"theta","minMu","rmax","initialization","mode",
			"multiclassPhiPoolingFunction","quiet","resultObject")
	
	
	if (resultObject!="S4"){
		res <- list(params, X.raw, X.all, norm$sizeFactors,
				AFinal, respFinal, fcFinal, INIFinal, pvalFinal,
				t(rINIT), t(meansINIT), pFinal, alphaFinal, 
				rFinal, meansFinal, sdsFinal)
		
		names(res) <- c("params","X","X.norm", "sizeFactors",
				"A","resp","logfc","INI","pval",
				"rINIT", "meansINIT","p","alpha",
				"r","means","sds")
		return(res)
	} else {
		res <- new("DEXUSResult")
		res@transcriptNames	 <- rownames(X.all)
		res@sampleNames		 <- colnames(X.all)
		res@inputData		 <- X.raw
		res@normalizedData   <- X.all  
		res@sizeFactors		 <- norm$sizeFactors
		res@INIValues        <- INIFinal
		res@INICalls		 <- INIFinal>0.1
		res@INIThreshold 	 <- 0.1
		res@pvals			 <- pvalFinal
		res@responsibilities <- respFinal	 
		res@posteriorProbs	 <- AFinal	
		res@logFC			 <- fcFinal 	
		res@conditionSizes	 <- alphaFinal	
		res@sizeParameters   <- rFinal		
		res@means            <- meansFinal			
		res@dispersions	     <- 1/rFinal		
		res@params           <- params 
		return(res)
	}
}
