#' get AUC for shared overlap nets
#' 
#' @details this predictor doesn't use GM to rank patient similarity. 
#' Rather a "called +" is a patient with CNVs that overlap feature-selected
#' nets
#' @param inFile (char) RR_changeNetSum_stats_den*Nets.txt. Contains a table
#' with pred/other OL for each score cutoff. Output of RR_featureTally()
#' @return (list) performance measures: 
#' 1) stats: (data.frame) ppv,f1,rec for each cutoff
#' 2) auc: (numeric) AUC-ROC
#' 3) prauc: (numeric) AUC-PR
#' @export
getAUC_OLnets <- function(inFile) {
	dat	<- read.delim(inFile,sep="\t",header=TRUE,as.is=TRUE)
	#dat[c(5,8,10:13)] <- dat[c(5,8,10:13)]/100

	comp_set <- data.frame(	
						score=dat$score,
						tp=dat$pred_ol,fp=dat$other_ol,
						# "-" that were correctly not called
						tn=dat$other_tot - dat$other_ol,
						# "+" that were not called 
						fn=dat$pred_tot - dat$pred_ol) 

	#plot(dat$other_pct,dat$pred_pct,
	#	  col='black', type="o",cex=0.5,lwd=2,
	#		xlim=c(0,1),ylim=c(0,1))

	stats <- netDx::perfCalc(comp_set)
	x <- stats$stats

	cat(sprintf("PRAUC = %1.2f; ",stats$prauc))
	cat(sprintf("ROCAUC = %1.2f\n", stats$auc))

	stats
}
