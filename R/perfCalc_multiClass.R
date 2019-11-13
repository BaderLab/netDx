#' Computes performance measures in multiway classification
#' 
#' @details First converts problem of N-way classification into N
#' binary classification problems ("A/not-A","B/not-B",etc.,).
#' Then computes precision/recall measure for each class. Finally returns
#' measures for each class as well as an overall average measure across
#' all binary.
#' @param realLabels (char) vector of real data labels
#' @param predLabels (char) vector of predicted data labels
#' @return (matrix) (N+1) named rows, one per class and one for average. 
#' Columns contain: tp,tn,fp,fn,ppv,recall,f1,accuracy
#' @export
perfCalc_multiClass <- function(realLabels,predLabels) {
	
	lvls <- unique(realLabels)
	out <- matrix(nrow=length(lvls)+1,ncol=8)
	colnames(out) <- c("tp","fp","tn","fn","ppv","recall","f1","acc")

	ctr <- 1
	for (lvl in lvls) {
		rLab <- realLabels; 
		rLab[setdiff(seq_len(length(rLab)),which(rLab%in%lvl))] <- "other"
		predLab <- predLabels; 
		predLab[setdiff(seq_len(length(rLab)),which(predLab%in%lvl))] <- "other"
		
		tp <- sum(rLab %in% lvl & rLab==predLab)
		tn <- sum(rLab %in% "other" & rLab==predLab)
		fp <- sum(rLab %in% "other" & rLab!=predLab)
		fn <- sum(rLab %in% lvl & rLab!=predLab)

		ppv	<- tp/(tp+fp)
		rec <- tp/(tp+fn)
		acc <- sum(rLab==predLab)/length(rLab)

		# F1 - harmonic mean of precision recall resolves to the 
		# formula below
		tp2 <- 2*tp
		f1 <- tp2/(tp2 + fp + fn)
		
		out[ctr,] <- c(tp,fp,tn,fn,ppv,rec,f1,acc)
		ctr <- ctr+1
	}
	out[nrow(out),] <- colMeans(out[-nrow(out),])
	
	out
}
