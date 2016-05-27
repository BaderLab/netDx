#' Nested cross validation for GM feature selection
#'
#' @details Runs two-tier cross validation by calling
#' \code{GM_runCV_featureSet()} for the inner loop. Computes patient
#' ranking and network tally
#' @param nOuterLoop (integer) number of outer loops for cross validation
#' @param setOloopSeed (integer) seed for random number generator that
#' in turn sets the seed for inner CV loop. The RNG samples numbers from
#' a vector 1:1000. Note: regardless, this value overrides the 
#' \code{setSeed} param of \code{makeCVqueries}.
#' @param outDir (char) directory to store results of GeneMANIA queries
#' @param pheno_DF (data.frame) patient ID and STATUS
#' @param predictClass (char) class (of pheno_DF$STATUS) being predicted
#' @param ... parameters for GM_runCV_featureSet()
#' @return (char) vector of paths to directories where result for each
#' outer loop is stored.
#' @export
GM_runCV_nested <- function(nOuterLoop=10L, setOloopSeed=42L,outDir,
	pheno_DF, predictClass,...) {
	if (!is.null(setOloopSeed)) {
		set.seed(setOloopSeed);
	}

	if(!file.exists(outDir)) dir.create(outDir)

	seedSet <- sample(1:1000, nOuterLoop,replace=FALSE)
	out <- list()
	for (seedLevel in seedSet) {
		fileSfx <- sprintf("CV_%i",seedLevel)
		newOut <- sprintf("%s/%s",outDir,fileSfx)
		if (!file.exists(newOut)) dir.create(newOut)

		GM_runCV_featureSet(setSeed=seedLevel,fileSfx=fileSfx,
							outDir=newOut,...)	

		# compute cross-validation error
		pFiles <- paste(newOut,dir(path=newOut,pattern="PRANK$"),sep="/")
		pFiles <- pFiles[grep(sprintf("^%s",fileSfx),basename(pFiles))]
		CV_perf(pFiles, pheno_DF, predictClass, newOut)

		# tally score for pathways
		nrankFiles <- paste(newOut,dir(path=newOut,pattern="NRANK$"),
							sep="/")
		pathwayRank	<- GM_networkTally(nrankFiles)
		write.table(pathwayRank,
					file=sprintf("%s/pathwayScore.txt",newOut),
					col=T,row=F,quote=F)

		out[[as.character(seedLevel)]] <- newOut
	}

	# compile results into a single table so we can take the average
	# across the outer loops
	# TODO change hard-coding of nrow
	tmp <- matrix(NA,nrow=10,ncol=nOuterLoop)
	for (k in 1:length(out)) {
		dat <- read.delim(sprintf("%s/CV_stats.txt",out[[k]]),h=T,as.is=T)
		tmp[,k] <- dat[,1]
	}
	AUC_CV <- tmp

	### now collect pathway score
	tmp <- list()
	for (k in 1:length(out)) {
		tmp[[k]] <- read.delim(sprintf("%s/pathwayScore.txt",out[[k]]),
						  h=T,as.is=T,sep=" ")
	}
	pscore <- tmp[[1]]
	for (k in 2:length(tmp)) {
		x <- merge(pscore,tmp[[k]],by='name',all.x=TRUE, all.y=TRUE,
			suffixes=c("",k));
		pscore <- x;
	}
	for (k in 2:ncol(pscore)) { 
			idx <- which(is.na(pscore[,k])); 
			pscore[idx,k] <- 0;
	}

	write.table(AUC_CV,file=sprintf("%s/AUC_nestedCV.txt", outDir),
				col=T,row=F,quote=F)
	write.table(pscore, 
				file=sprintf("%s/pathwayScore_nestedCV.txt",outDir),
				col=T,row=F,quote=F)

	return(list(AUC=AUC_CV, pscore=pscore))
}
