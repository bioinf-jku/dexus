# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>


#' @title Display the top ranked genes of the differential expression analysis.
#' 
#' @description Displays the result of DEXUS. A list of 
#' the 20 top-ranked genes is 
#' presented. In the unsupervised mode the genes are ranked by the I/NI call.
#' In case of two classes (groups/conditions) or multiple classes the 
#' genes are ranked by p-values.
#' @param object An instance of a "DEXUSResult".
#' @examples 
#' data(dexus)
#' r <- dexus(X[1:200, ])
#' show(r)
#' @return A data frame of the top ranked genes.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export
#' @noRd
#' @importFrom methods show

setMethod("show", "DEXUSResult",function(object){
			n <- min(10,nrow(object@normalizedData))
			#idx <- 1:n
			cat("Displaying the 10 top ranked genes of the analyis: \n")
			if (!is.null(object@params$mode)){
				if (object@params$mode=="unsupervised"){
					idx <- order(object@INIValues,decreasing=TRUE)[1:n]
					df <- data.frame(
							"Index"=idx,
							"Transcript"=object@transcriptNames[idx],
							"INIcall"=object@INICalls[idx],
							"INI"=round(object@INIValues[idx],3),
							"Mean"=round(object@means[idx, ,drop=FALSE],1))
					rownames(df) <- NULL
					print(df)
					
				} else  if (object@params$mode %in% c("two-class","multi-class")){
					idx <- order(object@pvals)[1:n]
					df <- data.frame(
							"Index"=idx,
							"Transcript"=object@transcriptNames[idx],
							"pvalues"=object@pvals[idx],
							"Mean"=round(object@means[idx, ,drop=FALSE],1))
					rownames(df) <- NULL
					print(df)
					
				} else {
					if (!all(is.na(object@pvals))){
						idx <- order(object@pvals)[1:n]
						df <- data.frame(
								"Index"=idx,
								"Transcript"=object@transcriptNames[idx],
								"pvalues"=object@pvals[idx],
								"Mean"=round(object@means[idx, ,drop=FALSE],1))
						rownames(df) <- NULL
						print(df)
						
					} else {
						if (!all(is.na(object@pvals))){
							idx <- order(object@INIValues,decreasing=TRUE)[1:n]
							df <- data.frame(
									"Index"=idx,
									"Transcript"=object@transcriptNames[idx],
									"INIcall"=object@INICalls[idx],
									"INI"=round(object@INIValues[idx],3),
									"Mean"=round(object@means[idx, ,drop=FALSE],1))
							rownames(df) <- NULL
							print(df)
							
						} else {
							message(paste("\"DEXUSResult\" contains neither",
											"p-values nor I/NI values."))
							
						}
					}
				}
			} else {
				if (!all(is.na(object@pvals))){
					idx <- order(object@pvals)[1:n]
					df <- data.frame(
							"Index"=idx,
							"Transcript"=object@transcriptNames[idx],
							"pvalues"=object@pvals[idx],
							"Mean"=round(object@means[idx, ,drop=FALSE],1))
					rownames(df) <- NULL
					print(df)
					
				} else {
					if (!all(is.na(object@INIValues))){
						idx <- order(object@INIValues,decreasing=TRUE)[1:n]
						df <- data.frame(
								"Index"=idx,
								"Transcript"=object@transcriptNames[idx],
								"INIcall"=object@INICalls[idx],
								"INI"=round(object@INIValues[idx],3),
								"Mean"=round(object@means[idx, ,drop=FALSE],1))
						rownames(df) <- NULL
						print(df)
					} else {
						message(paste("\"DEXUSResult\" contains neither",
										"p-values nor I/NI values."))
					}
				}
			}
			
			cat("\n")
			
			y <- length(object@INICalls)
			x <- length(which(object@INICalls))
			cat("Total number of transcripts: ",y,"\n")
			cat("Number of differentially expressed transcripts: ",x,"\n")
			cat("Percentage of differentially expressed transcripts: ",
					round(x/y,3)*100,"%","\n")
			
			#dexus::plot(object,n=min(20,nrow(object@normalizedData)))
		})


