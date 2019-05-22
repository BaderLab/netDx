#' Process GM PRANK files to get the ROC curve for the query
#'
#' @param pFile (char) path to PRANK file
#' @param pheno_DF (data.frame) patient IDs ("ID") and label("STATUS")
#' @param predClass (character) class label for which predictor is built
#' @param plotIt (logical) if TRUE plots ROC curve
#' @param verbose (logical) print messages
#' @export
#' @return (list) 
#' 1) predLbl: GeneMANIA scores (predicted labels). Higher score for
#' higher ranked patient. 
#' 2) realLbl: binary value indicating if patient label matches predictor
#' label (real labels)
#' 3) fullmat: pheno_DF merged with similarity scores ("similarityScore") and real label
#' ("isPredClass")
#' 4) roc: output of ROCRs performance(,"tpr","fpr") - ROC curve
#' 5) auc: output of ROCRs auc() 
#' 6) precall: output of ROCRs performance(, "prec","rec")
#' 7) f: output of ROCRs performance(,"f")
#' If < 2 patients in PRANK file, roc,auc, precall, f are all returned as
#' NA.
#' @examples
#' data(TCGA_mini)
#' prankFile <- sprintf("%s/extdata/GM_PRANK.txt", 
#' 	 path.package("netDx"))
#' x <- getPatientRankings(prankFile, pheno, "LumA")
getPatientRankings <- function(pFile,pheno_DF, predClass, plotIt=FALSE,
   verbose=FALSE) {
	# Read in PRANK file
	# need to skip comment line with new format
	dat <- read.table(pFile, sep="\t",header=TRUE, as.is=T, skip=1)

	pheno_DF$ID <- as.character(pheno_DF$ID)
	# 1 is what we predict, 0 is the other class
	pheno_DF$STATUS <- as.integer(pheno_DF$STATUS==predClass)

	if (verbose) cat(sprintf("%i total ; ", nrow(dat))) 
	dat	<- dat[which(!is.na(dat[,2])),]
	if (verbose) cat(sprintf("%i non-query entries in PRANK file\n", 
							 nrow(dat)))

	# match the pheno matrix to the labels
	midx <- match(dat[,1], pheno_DF$ID)
	if (all.equal(pheno_DF$ID[midx],dat[,1])!=TRUE) {
		cat("\t IDs in GM results don't match pheno\n"); browser()
	}
	curlbl <- pheno_DF[midx,,drop=FALSE]

	# ROCR prediction object requires the label assignments to range
	# from 0 to 1. GeneMANIA gene scores appear to be positive but 
	# unbounded.
	# We therefore rescale GM score to range from 0 to 1.
	curlbl <- cbind(curlbl, 
					similarityScore=order(dat[,2])/nrow(dat), 
					IsPredClass=curlbl$STATUS)
	curlbl <- curlbl[,-which(colnames(curlbl) %in% c("STATUS","TT_STATUS"))]

	# compute quality measures for classifier
	pred <- NA; perf <- NA; auc <- NA; precall <- NA; f <- NA
	# second condition is because ROCRs functions fail if there aren't
	# exactly two unique values in the real class label
	if (nrow(curlbl)>= 2 && length(unique(curlbl$IsPredClass))==2) {
		pred <- prediction(curlbl$similarityScore, curlbl$IsPredClass)
		perf <- performance(pred,"tpr","fpr")
		auc  <- performance(pred,"auc")@y.values[[1]]
		precall <- performance(pred,"prec","rec")
		f 		<- performance(pred,"f")
	
		if (plotIt) {
		plot(perf,main=sprintf("%i predictions; AUC= %1.2f",nrow(curlbl),auc),
			 bty='n',las=1,cex.axis=1.3)
		}
	}

	out <- merge(x=curlbl,y=pheno_DF,by="ID",all.y=TRUE)

	# return data incase caller would use it.
	return(list(predLbl=curlbl$similarityScore,realLbl=curlbl$IsPredClass,
				fullmat=out,pred=pred,
				roc=perf,auc=auc,
				precall=precall,f=f))
}
