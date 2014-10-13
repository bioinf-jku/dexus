\name{DEXUSResult-class}
\Rdversion{1.1}
\docType{class}
\alias{DEXUSResult-class}
\alias{conditionSizes,DEXUSResult-method}
\alias{dispersions,DEXUSResult-method}
\alias{INICalls,DEXUSResult-method}
\alias{INIThreshold,DEXUSResult-method}
\alias{INIThreshold<-,DEXUSResult-method}
\alias{INIValues,DEXUSResult-method}
\alias{inputData,DEXUSResult-method}
\alias{logFC,DEXUSResult-method}
\alias{means,DEXUSResult-method}
\alias{normalizedData,DEXUSResult-method}
\alias{params,DEXUSResult-method}
\alias{plot,DEXUSResult,missing-method}
\alias{posteriorProbs,DEXUSResult-method}
\alias{pvals,DEXUSResult-method}
\alias{responsibilities,DEXUSResult-method}
\alias{sampleNames,DEXUSResult-method}
\alias{show,DEXUSResult-method}
\alias{sizeFactors,DEXUSResult-method}
\alias{sizeParameters,DEXUSResult-method}
\alias{transcriptNames,DEXUSResult-method}
\alias{as.data.frame,DEXUSResult-method}
\alias{sort,DEXUSResult-method}
\alias{[,DEXUSResult,character-method}
\alias{[,DEXUSResult,numeric-method}
\alias{[,DEXUSResult,logical-method}
\alias{INI,DEXUSResult-method}


\title{Class \code{"DEXUSResult"}}
\description{
This class contains the result of an RNA-Seq data analysis. The class contains
the transcript names together with the parameters per condition, i.e., overdispersion 
and mean. Further it contains informative/non-informative values or p-values.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DEXUSResult", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{transcriptNames}}{The names of the transcripts, genes, exons, or regions of interest}
    \item{\code{sampleNames}}{The sample names as they were given in the input matrix.}
    \item{\code{inputData}}{The original read count matrix.}
    \item{\code{normalizedData}}{The normalized read count matrix.}
    \item{\code{sizeFactors}}{The size factors that were calculated for the normalization. This is that
factor that scales each column or sample.}
    \item{\code{INIValues}}{An informative/non-informative value for each sample that measures the evidence
for differential expression.}
    \item{\code{INIThreshold}}{The threshold for the I/NI values. Transcript with I/NI values above the
threshold will be considered as differentially expressed.}    
    \item{\code{INICalls}}{A binary value for each transcript indicating whether it is differentially expressed.}
    \item{\code{pvals}}{In case of two known conditions or multiple known conditions it is possible to
calculate a \eqn{p}-value for each transcript. This value is given in this slot.}
    \item{\code{responsibilities}}{A matrix of the size of the input matrix. It indicates the condition
for each sample and transcript. The condition named "1" is the major condition. All other
conditions are minor conditions. In case of supervised (two known conditions or multiple
known conditions) analyses this clustering matrix will be the same for all transcripts.}
    \item{\code{posteriorProbs}}{An array of the dimension of transcripts times samples times conditions.
It gives the probability that a certain read count x was generated under a condition.}
    \item{\code{logFC}}{The log foldchanges between the conditions. The reference is always condition "1".}
    \item{\code{conditionSizes}}{The ratio of samples belonging to that condition. 
    These are the \eqn{\alpha_i} values of the model.}
    \item{\code{sizeParameters}}{The size parameter estimates for each condition. These are the \eqn{r_i} values
of the model.}
    \item{\code{means}}{The mean of each condition. The \eqn{\mu_i} values of the model.}
    \item{\code{dispersions}}{The dispersion estimates for each condition. The inverse size parameters.}
    \item{\code{params}}{The input parameters of the DEXUS algorithm.}
  }
}
\section{Methods}{
  \describe{
    \item{[}{Subsetting of a DEXUSResult.}
    \item{as.data.frame}{Converts the result object into a data frame.}
    \item{conditionSizes}{Returns the condition sizes or 
    \eqn{\alpha_i} parameters of the model.}
    \item{dispersions}{Returns the dispersion, i.e. the inverse 
    size parameters, of the model. }
    \item{INI}{I/NI filtering of the result object.}
    \item{INICalls}{Returns a logical value 
    indication whether this transcript is differentially expressed or not.}
    \item{INIThreshold}{Returns the thresholds for the I/NI values.}
    \item{INIThreshold<-}{Sets the I/NI threshold. I/NI calls will be changed 
    accordingly.}
    \item{INIValues}{Returns the I/NI values.}
    \item{inputData}{Returns the input read counts.}
    \item{logFC}{Returns the log foldchange with respect to the first
    condition.}
    \item{means}{Returns the mean per condition.}
    \item{normalizedData}{Returns the normalized data.}
    \item{params}{Returns a list of input parameters of DEXUS.}
    \item{plot}{Plots a heatmap of the read counts of the top genes.}
    \item{posteriorProbs}{Returns an array of posterior probabilities.}
    \item{pvals}{Returns the \eqn{p}-values per transcript in supervised mode. }
    \item{responsibilities}{Returns the clustering vector.}
    \item{sampleNames}{Returns the sample names.}
    \item{show}{Displays a data frame of results.}
    \item{sizeFactors}{Returns the size factors used for normalization.}
    \item{sizeParameters}{Returns the size parameters, i.e. the \eqn{r_i} values
    of the model.}
    \item{sort}{Sorts the result object by I/NI values or p-values.}
    \item{transcriptNames}{Returns the transcript names.}
	 }
}
\author{
 Guenter Klambauer
}
\examples{
showClass("DEXUSResult")
}
\keyword{classes}
