# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>


#' @title Normalization of RNA-Seq count data.
#' 
#' @description Normalizes RNA-seq count data using previously published
#' approaches. Each samples'  read counts are corrected by a normalizing 
#' factor. The options are "RLE" by (Anders and Huber, 2010), 
#' and "upperquartile" by (Bullard et al., 2010).
#' 
#' @param X data a raw data matrix, where' columns are interpreted as samples 
#' and rows as genomic regions.
#' @param normalization method used for normalizing the reads. 
#' RLE is the method used by (Anders and Huber, 2010),
#' upperquartile is the Upper-Quartile method from (Bullard et al., 2010), and none
#' deactivates normalization. (Default = "RLE").
#' @return "list"  A list containing the normalized data (in its "X" component)
#' as well as the size-factors used for the normalization ("sizeFactors").
#' @references 
#' Anders, S. and Huber, W. (2010). \emph{Differential expression analysis
#' for sequence count data.} Genome Biol, 11(10), R106.
#' 
#' Bullard, J. H., Purdom, E., Hansen, K. D., and Dudoit, S. (2010).
#' \emph{Evaluation of statistical methods for normalization and 
#' differential expression in mRNA-seq experiments.} 
#' BMC Bioinformatics, 11, 94.
#' 
#' @examples 
#' data(dexus)
#' norm <- normalizeData(countsBottomly,"RLE")
#' 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas 
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export

normalizeData <- function(X, normalization){
	if (is.vector(X)){
		X <- matrix(X,nrow=1)
		X.int <- round(X)
		size.factors <- 1
	} else {
		X <- as.matrix(X)
		sf <- switch(normalization,
				RLE=normalize.rle(X),
				upperquartile=normalize.uq(X),
				none=rep(1, ncol(X)))
		X <- sweep(X, 2, sf, "/")
		X.int <- round(X)
	}
	mode(X) <- "numeric"
	return(list(X=X, sizeFactors=sf))
}

