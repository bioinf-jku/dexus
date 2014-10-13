# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>

# S4 class definition for the result object of a CNV detection method

#' @export
setClass("DEXUSResult",
		representation = representation
				(
				transcriptNames		= "character",
				sampleNames			= "character",
				inputData			= "matrix",
				normalizedData     	= "matrix",
				sizeFactors			= "numeric",
				INIValues        	= "numeric",
				INIThreshold		= "numeric",
				INICalls			= "logical",
				pvals				= "numeric",
				responsibilities 	= "matrix",
				posteriorProbs		= "array",
				logFC				= "matrix",
				conditionSizes		= "matrix",
				sizeParameters 		= "matrix",
				means				= "matrix",
				dispersions			= "matrix",
				params				= "list"),
		prototype = prototype
		(
				transcriptNames		= character(),
				sampleNames			= character(),
				inputData			= matrix(),
				normalizedData     	= matrix(),
				sizeFactors			= numeric(),
				INIValues        	= numeric(),
				INIThreshold 		= 0.1,
				INICalls			= logical(),
				pvals 				= numeric(),
				responsibilities 	= matrix(),
				posteriorProbs		= array(),
				logFC				= matrix(),
				conditionSizes		= matrix(),
				sizeParameters 		= matrix(),
				means				= matrix(),
				dispersions			= matrix(),
				params				= list()
				)
)
