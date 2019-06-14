#' Pick best cutoff following resampling of training data
#'
#' @param perf_resamp (list) the output of evalPerfOverCutoffs(), 
#' with one entry for each round of resampling
#' @return (list) with performance at each cutoff:
#' 1) perfmat (matrix): tpr, fpr, accuracy,ppv for each resampling. One 
#' row per cutoff, 4 columns per round of resampling.
#' 2) bestCutoff (integer): cutoff with the highest accuracy
#' 3) mu (matrix): average of tpr, fpr, accuracy, ppv across resamplings.
#' one row per cutoff.
#' @export
resampling_pickBestCutoff <- function(perf_resamp) {

cat("* Picking best cutoff from resampling test\n")

perfmat <- matrix(NA, 
	nrow=nrow(perf_resamp[[1]][["confmat"]]),
	ncol=4*length(perf_resamp))
colnames(perfmat) <- rep(c("tpr","fpr","accuracy","ppv"),
	 length(perf_resamp))

ctr <- 1
for (k in 1:length(perf_resamp)) {
	perfmat[,ctr:(ctr+3)] <- perf_resamp[[k]][["confmat"]][,6:9]
	ctr <- ctr+4
}

mu <- matrix(NA,nrow=nrow(perfmat),ncol=4)
for (k in 1:4) {
	print(seq(k,ncol(perfmat),4))
	mu[,k] <- rowMeans(perfmat[,seq(k,ncol(perfmat),4)])
}
colnames(mu) <- c("tpr","fpr","accuracy","ppv")

# cutoff is the one with best accuracy
bestIdx <- which.max(mu[,3]) 
bestCutoff <- perf_resamp[[1]][["confmat"]][bestIdx,1]
cat(sprintf("Best cutoff = %i ; TPR=%1.2f; FPR=%1.2f; Accuracy=%1.2f; PPV=%1.2f\n",
	bestCutoff, mu[bestIdx,1],mu[bestIdx,2],mu[bestIdx,3],mu[bestIdx,4]))

return(list(perfmat=perfmat, bestCutoff=bestCutoff, mu=mu))

}
