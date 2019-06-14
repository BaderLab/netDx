#' pick best cutoff for CNV-based predictor
#'
#' @param resampPerf (matrix) output of RR_featureTally stored in 
#' <outDir>/resampPerf.Rdata. Has one row per score, and (7*resamp) columns. 
#' 7 columns for each resampling round. see RR_featureTally() for more details.
#' scores must be in rownames(resampPerf)
#' @return (list)
#' 1) bestCutoff (integer) score with highest accuracy
#' 3) accuracy (matrix) accuracy for each score cutoff
#' @export
resampling_pickBestCutoff_CNV <- function(resampPerf) {

	numresamp <- ncol(resampPerf)/7
	resampPerf <- as.data.frame(resampPerf)
	
	acc <- matrix(nrow=nrow(resampPerf),ncol=numresamp)
	rownames(acc) <- rownames(resampPerf)
	for (k in 1:numresamp) {
		tpc		<- sprintf("pred_OL_%i",k)
		posc	<- sprintf("pred_total_%i",k) 
		fpc		<- sprintf("other_OL_%i",k)
		negc 	<- sprintf("other_total_%i",k)

		# tn = (fp+tn)-fp
		tn <- resampPerf[,negc]-resampPerf[,fpc]		

		num <- resampPerf[,tpc] + tn
		den <- resampPerf[,posc] + resampPerf[,negc]
		acc[,k] <- num/den
	}

	num_x <- nrow(resampPerf)
	plot(0,0,xlim=c(0,num_x+1),ylim=c(0,1),type='n',
		xlab="score",ylab="accuracy",
		main="CNV predictor\nAccuracy over resamplings")
	clrs <- brewer.pal(name="Dark2",n=ncol(acc))
	for (k in 1:ncol(acc)) {
		points(1:num_x,acc[,k],col=clrs[k],
			type='o',pch=16)
	}
		points(1:num_x,rowMeans(acc),col="black",
			type='o',pch=4,lwd=2)
	legend("topleft",legend=c(1:numresamp,'mean'),
		fill=NA,border=NA,lty=1,lwd=1,col=c(clrs,'black'),
		pch=c(rep(16,length(clrs)),4),bty='n')
	
	acc <- cbind(acc,rowMeans(acc))
	idx <- which.max(acc[,4])
	bestCutoff <- as.integer(rownames(acc)[idx])
	cat(sprintf("Best mean resampling accuracy at score = %i; acc= %1.1f%%\n",
		bestCutoff, acc[idx,4]*100))
	colnames(acc) <- c(paste("Resampling", 1:numresamp,sep="_"),"Mean")
	
	return(list(accuracy=acc, bestCutoff=which.max(acc[,4])))
}
