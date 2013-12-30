# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>

#' @title Accessors for a "DEXUSResult".
#' @name accessors
#' 
#' @aliases transcriptNames sampleNames inputData normalizedData 
#' sizeFactors INIValues INIThreshold INICalls pvals responsibilities
#' posteriorProbs logFC conditionSizes sizeParameters means dispersions
#' params
#' 
#' @rdname accessors
#' 
#' @description These generic functions return the slots
#' of an RNA-Seq analysis performed by DEXUS. 
#' The results of DEXUS are stored as an instance of 
#' \code{\link{DEXUSResult-class}}.
#' 
#' @param object An instance of "DEXUSResult".
#' @examples
#' data(dexus)
#' result <- dexus(countsBottomly[1:20,1:10])
#' transcriptNames(result)
#' sampleNames(result)
#' inputData(result)
#' normalizedData(result)
#' sizeFactors(result)
#' INIValues(result)
#' INIThreshold(result)
#' INICalls(result)
#' pvals(result)
#' responsibilities(result)
#' posteriorProbs(result)
#' logFC(result)
#' conditionSizes(result)
#' sizeParameters(result)
#' means(result)
#' dispersions(result)
#' params(result)
#' 
#' @export transcriptNames 
#' @export sampleNames 
#' @export inputData 
#' @export normalizedData 
#' @export sizeFactors
#' @export INIValues 
#' @export INIThreshold 
#' @export INICalls 
#' @export pvals 
#' @export responsibilities
#' @export posteriorProbs 
#' @export logFC 
#' @export conditionSizes
#' @export sizeParameters 
#' @export means 
#' @export dispersions
#' @export params
#' @return The accessor functions return a the matrices or vectors contained
#' in the corresponding slot of the "DEXUSResult".
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}

setMethod("transcriptNames", signature = "DEXUSResult", 
		definition = function(object) object@transcriptNames)

setMethod("sampleNames", signature = "DEXUSResult", 
		definition = function(object) object@sampleNames)

setMethod("inputData", signature = "DEXUSResult", 
		definition = function(object) object@inputData)

setMethod("normalizedData", signature = "DEXUSResult", 
		definition = function(object) object@normalizedData)

setMethod("sizeFactors", signature = "DEXUSResult", 
		definition = function(object) object@sizeFactors)

setMethod("INIValues", signature = "DEXUSResult", 
		definition = function(object) object@INIValues)


setMethod("INIThreshold", signature = "DEXUSResult", 
		definition = function(object) object@INIThreshold)

setMethod("INICalls", signature = "DEXUSResult", 
		definition = function(object) object@INICalls)

setMethod("pvals", signature = "DEXUSResult", 
		definition = function(object) object@pvals)

setMethod("responsibilities", signature = "DEXUSResult", 
		definition = function(object) object@responsibilities)

setMethod("posteriorProbs", signature = "DEXUSResult", 
		definition = function(object) object@posteriorProbs)

setMethod("logFC", signature = "DEXUSResult", 
		definition = function(object) object@logFC)


setMethod("conditionSizes", signature = "DEXUSResult", 
		definition = function(object) object@conditionSizes)

setMethod("sizeParameters", signature = "DEXUSResult", 
		definition = function(object) object@sizeParameters)

setMethod("means", signature = "DEXUSResult", 
		definition = function(object) object@means)

setMethod("dispersions", signature = "DEXUSResult", 
		definition = function(object) object@dispersions)

setMethod("params", signature = "DEXUSResult", 
		definition = function(object) object@params)



#' @title Set the I/NI threshold.
#' @aliases INIThreshold<-
#' @aliases INIThreshold-set

#' 
#' @description This generic function sets the 
#' threshold of the I/NI value. Transcripts with I/NI values above the
#' I/NI threshold are considered as differentially expressed. 
#' The results of DEXUS are stored as an instance of 
#' \code{\link{DEXUSResult-class}}.
#' 
#' @param object An instance of "DEXUSResult".
#' @param value A numeric to be used for thresholding the I/NI values.
#' @examples
#' data(dexus)
#' result <- dexus(countsBottomly[1:20,1:10])
#' INIThreshold(result) <- 0.1
#' @return \code{INIThreshold<-} returns an instance of "DEXUSResult".
#' @name INIThreshold<-
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export 
setReplaceMethod("INIThreshold", signature="DEXUSResult",
		definition=function(object,value)
		{
			slot(object,"INIThreshold") <- value
			slot(object,"INICalls") <- object@INIValues > value
			y <- length(object@INICalls)
			x <- length(which(object@INICalls))
			cat("Total number of transcripts: ",y,"\n")
			cat("Number of differentially expressed transcripts: ",x,"\n")
			cat("Percentage of differentially expressed transcripts: ",
					round(x/y,3)*100,"%","\n")
			
			object
		}
)



#' @title Obtaining the result data frame.
#' 
#' @description This generic function extracts the result of the analysis
#' from the "DEXUSResult" class as a data frame.
#' 
#' The results of DEXUS are stored as an instance of 
#' \code{\link{DEXUSResult-class}}.
#' 
#' @param object An instance of "DEXUSResult".
#' @examples
#' data(dexus)
#' result <- dexus(countsBottomly[1:20,1:10])
#' as.data.frame(result)
#' @return \code{as.data.frame} returns a data frame.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @noRd
#' @exportMethod as.data.frame


setMethod("as.data.frame", signature = c(x="DEXUSResult"), 
		definition = function(x){
			
			n <- nrow(x@normalizedData)
			idx <- 1:n
			df <- data.frame(
					"Transcript"=x@transcriptNames[idx],
					"INIcall"=x@INICalls[idx],
					"INI"=x@INIValues[idx],
					"pval"=x@pvals[idx],
					"Mean"=x@means[idx, ,drop=FALSE],
					"logFC"=x@logFC[idx, ,drop=FALSE],
					"conditionSize"=x@conditionSizes[idx, ,drop=FALSE],
					"dispersion"=x@dispersions[idx, ,drop=FALSE],
					"responsibilities"=x@responsibilities[idx, ,drop=FALSE]
			)
			rownames(df) <- NULL
			return(df)
		})



#' @title Subsetting a "DEXUSResult".
#' 
#' @description Information about specific transcripts can be accessed in 
#' the "DEXUSResult" object by using the standard brackets "[idx]" for
#' subsetting. Either transcript names or transcript indices can be used.
#' 
#' @name `[`
#' @aliases `[`,DEXUSResult,numeric-method
#' @aliases `[`,DEXUSResult,character-method
#' @aliases `[`,DEXUSResult,logical-method

#' @param x "DEXUSResult"
#' @param i Either a numeric vector of indices or a character vector containing
#' the transcript names.
#' @return An instance of "DEXUSResult".
#' 
#' @examples
#' data(dexus)
#' res <- dexus(countsBottomly[1:100, ])
#' res["ENSMUSG00000000486"]
#' res[50:55]
#' 
#' @rdname DEXUSResult-subset
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export


setMethod("[", signature = c(x="DEXUSResult",i="numeric"), 
		definition = function(x,i){
			
			n <- length(x@transcriptNames)
			y <- new("DEXUSResult")
			
			if (!all(i %in% 1:n)){
				stop(paste("Given Indices outside the range of the",
								"transcripts' indices."))
			}  else {
				y@transcriptNames	= x@transcriptNames[i]
				y@sampleNames		= x@sampleNames
				y@inputData			= x@inputData[i, ,drop=FALSE]
				y@normalizedData    = x@normalizedData[i, ,drop=FALSE]
				y@sizeFactors		= x@sizeFactors
				y@INIValues        	= x@INIValues[i]
				y@INIThreshold 		= x@INIThreshold
				y@INICalls			= x@INICalls[i]
				y@pvals 			= x@pvals[i]
				y@responsibilities 	= x@responsibilities[i, ,drop=FALSE]
				y@posteriorProbs	= x@posteriorProbs[, ,i,drop=FALSE]
				y@logFC				= x@logFC[i, ,drop=FALSE]
				y@conditionSizes	= x@conditionSizes[i, ,drop=FALSE]
				y@sizeParameters 	= x@sizeParameters[i, ,drop=FALSE]
				y@means				= x@means[i, ,drop=FALSE]
				y@dispersions		= x@dispersions[i, ,drop=FALSE]
				y@params			= x@params
			}
			
			return(y)
		})

setMethod("[", signature = c(x="DEXUSResult",i="logical"), 
		definition = function(x,i){
			
			n <- length(x@transcriptNames)
			y <- new("DEXUSResult")
			if (length(i)!=n){
				stop(paste("Given logical vector must have as many elements",
								"as transcripts in the result object."))
			}
			
			i <- which(i)		
			
			y@transcriptNames	= x@transcriptNames[i]
			y@sampleNames		= x@sampleNames
			y@inputData			= x@inputData[i, ,drop=FALSE]
			y@normalizedData    = x@normalizedData[i, ,drop=FALSE]
			y@sizeFactors		= x@sizeFactors
			y@INIValues        	= x@INIValues[i]
			y@INIThreshold 		= x@INIThreshold
			y@INICalls			= x@INICalls[i]
			y@pvals 			= x@pvals[i]
			y@responsibilities 	= x@responsibilities[i, ,drop=FALSE]
			y@posteriorProbs	= x@posteriorProbs[, ,i,drop=FALSE]
			y@logFC				= x@logFC[i, ,drop=FALSE]
			y@conditionSizes	= x@conditionSizes[i, ,drop=FALSE]
			y@sizeParameters 	= x@sizeParameters[i, ,drop=FALSE]
			y@means				= x@means[i, ,drop=FALSE]
			y@dispersions		= x@dispersions[i, ,drop=FALSE]
			y@params			= x@params
			
			
			return(y)
		})


setMethod("[", signature = c(x="DEXUSResult",i="character"), 
		definition = function(x,i){
			
			y <- new("DEXUSResult")
			
			if (!all(i %in% x@transcriptNames)){
				stop("Given transcript name not in result object.")
			}  else {
				ii <- match(i,x@transcriptNames)
				y@transcriptNames	= x@transcriptNames[ii]
				y@sampleNames		= x@sampleNames
				y@inputData			= x@inputData[i, ,drop=FALSE]
				y@normalizedData    = x@normalizedData[i, ,drop=FALSE]
				y@sizeFactors		= x@sizeFactors
				y@INIValues        	= x@INIValues[i]
				y@INIThreshold 		= x@INIThreshold
				y@INICalls			= x@INICalls[i]
				y@pvals 			= x@pvals[i]
				y@responsibilities 	= x@responsibilities[i, ,drop=FALSE]
				y@posteriorProbs	= x@posteriorProbs[, ,ii,drop=FALSE]
				y@logFC				= x@logFC[i, ,drop=FALSE]
				y@conditionSizes	= x@conditionSizes[i, ,drop=FALSE]
				y@sizeParameters 	= x@sizeParameters[i, ,drop=FALSE]
				y@means				= x@means[i, ,drop=FALSE]
				y@dispersions		= x@dispersions[i, ,drop=FALSE]
				y@params			= x@params
			}
			
			return(y)
		})




#' @title Sorting a DEXUS result.
#' @description This function sorts the result object by I/NI values or
#' p-values such that the transcripts with the highest I/NI value or the lowest
#' p-value are ranked first.
#' 
#' @name sort
#' 
#' @param object An instance of "DEXUSResult".
#' @return An instance of "DEXUSResult".
#' 
#' @examples
#' data(dexus)
#' res <- dexus(countsBottomly[1:100, ])
#' sort(res)
#' 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @exportMethod sort

setMethod("sort", 
		signature= c(x="DEXUSResult"),
		definition=function(x){
			n <- nrow(x@normalizedData)
			
			if (!is.null(x@params$mode)){
				if (x@params$mode=="unsupervised" | x@params$mode=="semi-supervised"){
					idx <- order(x@INIValues,decreasing=TRUE)[1:n]		
				} else  if (x@params$mode %in% c("two-class","multi-class")){
					idx <- order(x@pvals)[1:n]	
				} else {
					if (!all(is.na(x@pvals))){
						idx <- order(x@pvals)[1:n]						
					} else {
						if (!all(is.na(x@pvals))){
							idx <- order(x@INIValues,decreasing=TRUE)[1:n]
						} else {
							stop(paste("\"DEXUSResult\" contains neither",
											"p-values nor I/NI values."))							
						}
					}
				}
			} else {
				if (!all(is.na(x@pvals))){
					idx <- order(x@pvals)[1:n]
				} else {
					if (!all(is.na(x@INIValues))){
						idx <- order(x@INIValues,decreasing=TRUE)[1:n]						
					} else {
						stop(paste("\"DEXUSResult\" contains neither",
										"p-values nor I/NI values."))
					}
				}
			}

			return(x[idx])
		})



#' @title I/NI filtering of a DEXUS result.
#' @description This function filters the result object for informative 
#' transcripts. Transcripts with an I/NI value below the given threshold
#' are filtered out.
#' 
#' @name INI
#' 
#' @param object An instance of "DEXUSResult".
#' @param threshold A numeric determining the threshold for the I/NI values.
#' @return An instance of "DEXUSResult".
#' 
#' @examples
#' data(dexus)
#' res <- dexus(countsBottomly[1:100, ])
#' INI(res)
#' 
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and 
#' Thomas Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export



setMethod("INI", signature="DEXUSResult",
		definition=function(object,threshold=0.1)
		{
			if (all(is.na(object@INIValues)))
				stop(paste("No I/NI values found in the result object.",
								"INI filter can only be applied if DEXUS",
								"was run in unsupervised mode."))
			
			idx <- object@INIValues > threshold
			y <- length(idx)
			x <- length(which(idx))
			cat("Total number of transcripts: ",y,"\n")
			cat("Number of differentially expressed transcripts: ",x,"\n")
			cat("Percentage of differentially expressed transcripts: ",
					round(x/y,3)*100,"%","\n")
			cat("\n")
			
			return(sort(object[which(idx)]))
		}
)
