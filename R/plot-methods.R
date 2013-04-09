# Copyright (C) 2013 Guenter Klambauer and Thomas Unterthiner
# <klambauer@bioinf.jku.at>

#' @title Visualization of a result of the DEXUS algorithm.
#' 
#' @description Plots a heatmap of the log read counts of the top ranked genes
#' or of selected genes.
#' @param x An instance of "CNVDetectionResult" 
#' @param idx The indices or the transcript names of the transcripts that 
#' should be visualized as heatmap.
#' @param cexSamples Size of the column labels, i.e. the samples.
#' @param cexGenes Size of the row labels, i.e. the transcripts.
#' @param newColNames renames the samples.
#' @param type Mark the samples, that do not belong to the major class by
#' crosses ("crosses"), or boxes ("boxes"). 
#' @examples
#' data(dexus)
#' r <- dexus(countsBottomly[1:100, ])
#' plot(r)
#' @return Generates a heatmap of the expression values of the top-ranked transcripts.
#' @author Guenter Klambauer \email{klambauer@@bioinf.jku.at} and Thomas
#' Unterthiner \email{unterthiner@@bioinf.jku.at}
#' @export
#' @importFrom graphics plot

setMethod("plot", signature(x="DEXUSResult",y="missing"),
		function(x, idx, cexSamples=0.5,cexGenes=1,
				newColNames=NULL,type="crosses",cexCrosses=2){
			library(RColorBrewer)			
			nmax <- length(x@INIValues)
			if (missing(idx))
				idx <- order(x@INIValues,decreasing=TRUE)[1:min(15,nmax)]
			
			idx <- rev(idx)
			xf <- x[idx]
			
			
			X <- log(xf@normalizedData+0.01)
			if (!is.null(newColNames)){
				colnames(X) <- newColNames
			} 
			rownames(X) <- xf@transcriptNames
			
			m <- ncol(X)
			n <- nrow(X)
			resp <- t(xf@responsibilities)
			
			X <- t(X)
			
			par(oma=c(1,6,0,0))
			rgb.palette <- colorRampPalette(c("white", "darkblue"),space = "rgb")
			image(X,xaxt="n",yaxt="n",col = rgb.palette(100))
			if (n==1){
				ia <- 1
				yy <- c(-1,1)
			} else {
				ia <- 1/(n-1)
				yy <- seq(-ia/2,1+ia/2,ia)
			}
			
			ib <- 1/(m-1)
			xx <- seq(-ib/2,1+ib/2,ib)
			for (j in 1:n){
				for (i in 1:m){
					if (type=="innerboxes"){
						if (resp[i,j]!=1){
							rect(xleft=xx[i]+ib/15,xright=xx[i+1]-ib/15,
									ybottom=yy[j]+ia/15,ytop=yy[j+1]-ia/15,
									lwd=5,border="red")
						}
					} else if (type=="boxes"){
						if (resp[i,j]!=1){
							rect(xleft=xx[i],xright=xx[i+1],
									ybottom=yy[j],ytop=yy[j+1],
									lwd=2,border="red")
						}
					} else {
						if (resp[i,j]!=1){
							points(x=(xx[i]+xx[i+1])/2,y=(yy[j]+yy[j+1])/2,
									pch=4+(resp[i,j]-2),col="red",
									lwd=3,cex=cexCrosses)
						}
					}
				}	
			}
			axis(1,at=seq(0,1,ib),labels=rownames(X),las=2,cex.axis=cexSamples)
			if (n==1){
				axis(2,at=0,labels=colnames(X),las=2,cex.axis=cexGenes)
			} else {
				axis(2,at=seq(0,1,ia),labels=colnames(X),las=2,cex.axis=cexGenes)
			}
		})

