#' Computes variety of predictor evaluation measures based on the confusion
#' matrix
#'
#' @param dat (data.frame): 5 columns: score, tp, fp, tn, fn. 
#' One row per cutoff
#' score for feature selection
#' @return (list)
#' stats (data.frame): score, f1, ppv, precision and recall. One row
#' per cutoff for feature selection
#' auc (numeric between 0 and 1): AUC of overall ROC curve
#' prauc (numeric between 0 and 1): AUC of overall precision-recall curve
#' @importFrom pracma trapz
#' @examples
#' data(confmat)
#' x <- perfCalc(confmat)
#' @export
perfCalc <- function(dat) {
	dat <- na.omit(dat)
	# F1 - harmonic mean of precision recall resolves to the formula below
	tp2 <- 2*dat$tp
	f1 <- tp2/(tp2 + dat$fp + dat$fn)

	#precision recall curve

	# precision = positive predictive value (pr = ppv)
	ppv	<- dat$tp/(dat$tp+dat$fp)
	rec <- dat$tp/(dat$tp+dat$fn)
	# trapz integrates from right to left, so you need to apply rev()
	# otherwise you get a negative area.
	prauc <- pracma::trapz(rev(rec),rev(ppv))

	#roc auc
	x <- dat$fp/(dat$fp+dat$tn)
	y <- dat$tp/(dat$tp+dat$fn)
	
	x <- c(0,rev(x),1); y <- c(0,rev(y),1)
	auc <- pracma::trapz(x,y)
	out <- data.frame(score=dat$score,ppv=ppv,f1=f1,rec=rec)

	return(list(stats=out,auc=auc,prauc=prauc))
}
