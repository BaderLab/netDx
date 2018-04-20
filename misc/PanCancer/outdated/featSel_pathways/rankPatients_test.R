#' snippet of code from GM_getQueryROC 
#' @param pFile PRANK file
#' @param pheno_DF (data.frame) must have ID, STATUS,TT_STATUS column
rankPatient <- function(pFile,pheno_DF,verbose=FALSE) {
 	dat <- read.table(pFile, sep="\t",header=TRUE, as.is=T)
	if (verbose) cat(sprintf("%i total ; ", nrow(dat)))
  dat <- dat[which(!is.na(dat[,2])),]
  if (verbose) cat(sprintf("%i non-query entries in PRANK file\n",
               nrow(dat)))

  # ROCR prediction object requires the label assignments to range
  # from 0 to 1. GeneMANIA gene scores appear to be positive but
  # unbounded.
  # We therefore rescale GM score to range from 0 to 1.
  dat <- cbind(dat,
          GM_score=order(dat[,2])/nrow(dat),
					GM_weight=dat[,2])
	colnames(dat)[1] <- "ID"
	return(dat)
}
