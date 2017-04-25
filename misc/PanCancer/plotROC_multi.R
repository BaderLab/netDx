#' plot multiple ROC curves with average.
#'
#' Plots average ROC curve with individual ROC curves imposed.
#' @param (list) ROCR::performance objects
#' @param (numeric) index to highlight in a bold colour
plotROC_multi <- function(inList,idx2show) { 
	par(las=1,cex.axis=1.3,cex.lab=1.3)
	plot(0,0,type='n',xlim=c(0,1),ylim=c(0,1),xlab="FPR",ylab="TPR",
		bty='n')
	for (k in 1:length(inList)) {
		x <- inList[[k]]@x.values[[1]]
		y <- inList[[k]]@y.values[[1]]
		if (k==idx2show)
			points(x,y,type="l",col="red",lwd=2)
		else 
			points(x,y,type="l",col=rgb(0.5,0.5,0.5,0.3))

	}	
	abline(0,1,col='darkblue',lwd=2)
}
