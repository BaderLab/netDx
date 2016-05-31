#' ROC curve for cross validation runs
#'
#' @param fList (char) list of GM result files
#' @param pheno_DF (data.frame) ID: patient ID ; STATUS: patient STATUS
#' @param predClass (char) class for which predictor is being built
#' @param outDir (char) path to output files
#' @param justGetError (logical)  if FALSE , computes all stats (plots AUC 
#' curve etc., etc., if TRUE, simply returns cross-validation error.
#' The latter is useful for iterative feature selection.
#' @param errorType (char) one of [AUPR|maxF1]; the type of performance 
#'	measure to use for evaluating cross-validation
#' @param verbose (logical) print messages
#' @return Depends on value of \code{justGetError}. If justGetError=TRUE,
#' then returns the statistic requested by the \code{errorType} parameter.
#' Note: if any of the folds returns an NA the prediction/roc/auc/f, which
#' happens if only one datapoint is available in the PRANK file, then 
#' the error metric returned is 0. 
#' Else returns nothing, but has the side effect of printing the AUC
#' for all folds to the output file CV_stats.txt in \code{outDir}. In 
#' this scenario also plots the n-fold AUC curve in \code{CV_AUC.pdf}.
#' @export
#' @import ROCR
#' @importFrom pracma trapz
CV_perf <- function(fList,pheno_DF,predClass,outDir,justGetError=FALSE,
	errorType="AUPR",verbose=TRUE) {
if (! errorType %in% c("AUPR","maxF1")) 
		stop("errorType must be AUPR or maxF1")

predCV	<- list()
lblCV	<- list()
pred    <- list()
perf    <- list()
auc     <- list()
fstat   <- list()
precall <- list()


ctr		<- 1;
if (!justGetError) {
	pdf(sprintf("%s/GM_patientRanking.pdf",outDir),width=8,height=8)
	par(mfrow=c(2,2),las=1,font=2,cex.axis=1.3,bty='n')
}

#scoreList <- list() # GM patient scores
for (curF in fList) {
	if (verbose) cat(sprintf("%s => ", basename(curF)))
	
	dat <- read.table(curF, sep="\t",header=T,as.is=T)
	dat	<- dat[which(!is.na(dat[,2])),]
	res <- GM_getQueryROC(curF, pheno_DF, predClass,verbose=verbose)
	
	# ROCR prediction object requires the label assignments to range
	# from 0 to 1. GeneMANIA gene scores appear to be positive but 
	# unbounded.
	# We therefore rescale GM score to range from 0 to 1.
	predCV[[ctr]]	<- res$predLbl
	lblCV[[ctr]]	<- res$realLbl
	pred[[ctr]] 	<- res$pred
	perf[[ctr]]		<- res$roc
	auc[[ctr]]		<- res$auc
	precall[[ctr]]	<- res$precall
	fstat[[ctr]]	<- res$f

	if (length(res$predLbl)==1) {
		if (!is.na(res$predLbl)) {
			if (!justGetError) {
			clrs <- rep(rgb(0,0,0,0.3),nrow(dat))
			clrs[which(res$isPredClass>0)] <- "red"
			plot(1:length(dat[,2]),dat$Score,
				 col=clrs,pch=16,cex=0.6,
				 xlab="PRANK", ylab="Patient Score",
				 main=basename(curF),cex.main=0.7)
			legend("topright",fill=NA,border=NA,col=c("black","red"),
				   pch=16,legend=c("not in class", "in class"),bty='n')
			}
	}}
	ctr	<- ctr+1
}
if (!justGetError) dev.off()

######
###### now compute model performance measures

#pred	<- prediction(predCV, lblCV)
#perf	<- performance(pred,"tpr","fpr")
#tmp		<- performance(pred, "auc"); 
#f		<- performance(pred,"f")
#precall	<- performance(pred,"prec","rec")
#

# just return AUPR
if (justGetError) {
	if (errorType %in% "AUPR") myarr <- precall else myarr <- f
	ln 	<- length(myarr)
	tmp	<-	sapply(1:ln, function(i) {
		out <- 0
		isEmpty <- suppressWarnings(is.na(myarr[[i]]))
		if (!isEmpty){
			if (errorType %in% "AUPR") {
				out <- trapz(myarr[[i]]@x.values[[1]][-1],
						 myarr[[i]]@y.values[[1]][-1])
			} else {
				maxInd <- which.max(myarr[[i]]@y.values[[1]])
				out		<- myarr[[i]]@y.values[[1]][maxInd]
			}
		}
		out
	})
	if (verbose) print(tmp)
	return(mean(tmp))
}

# auc
cat("AUC\n----------------------\n")
auc 	<- unlist(auc)
print(auc)

out_stats	<- matrix(NA,nrow=length(auc),ncol=4)
colnames(out_stats) <- c("AUC","max-F1","cutoff for max-F1","AUPR")
out_stats[,1]	<- auc

# precision recall
cat("Precision recall F1\n----------------------\n")
for (i in 1:nrow(out_stats)) {
	maxInd	<- which.max(fstat[[i]]@y.values[[1]])
	fmax	<- fstat[[i]]@y.values[[1]][maxInd]
	cutoff	<- fstat[[i]]@x.values[[1]][maxInd]
	# NOTE while computing AUPR. The first y values in the PR curve
	# returned by ROCR is always NAN. If we include this, the AUPR
	# is always NAN. So in the trapz() call below we exclude this
	# first element.
	out_stats[i,2:4] <- c(fmax,cutoff,
						 trapz(precall[[i]]@x.values[[1]][-1],
							   precall[[i]]@y.values[[1]][-1])
						 )
}
cat("\nmax F1 for each fold:\n")
print(out_stats[,2])

write.table(out_stats, file=sprintf("%s/CV_stats.txt",outDir), 
		sep="\t",col=T,row=F,quote=F)

}
