#' Build predictor with n-way resampling
#'
#' @details Builds a binary classifier, with feature selection that 
#' resamples training data. Current limitations. Binary prediction only. 
#' Patient nets must
#' be limited to full matrix data. Future extensions should allow some or 
#' all RangeSet nets. NOTE: To ensure reproducibility, set the RNG seed
#' before running this method. This method does not currently allow users
#' to explicitly set the RNG seed for specific functions 
#' (e.g. splitTestTrain)
#' @param pheno (data.frame) patient IDs (ID) and class of interest (STATUS)
#' @param pdat (matrix) patient data. Rows are patients, columns are unit
#' measures (e.g. genes)
#' @param predClass (char) class to build predictor for. Must be a value
#' in the STATUS column of the pheno param
#' @param unitSets (list) unit groupings, each of which will be converted to
#' its own patient similarity network. 
#' @param pctT (numeric 0 to 1) fraction of samples for training
#' @param numResamples (integer) number of resamplings for cross validation
#' @param nFoldCV (integer) number of folds for cross-validation within each
#' resampling
#' @param numCores (integer) number of parallel processing jobs
#' @param GMmemory (integer) memory (Gb) for GeneMANIA to run.
#' @param outDir (char) path to output directory
#' @param overwrite (logical) Overwrite existing results in outDir? By 
#' default, no, return warning; but set to TRUE if outDir should be deleted
#' and populated with fresh results from this function.
#' @param ... parameters for makePSN_NamedMatrix()
#' @return (list) Includes:
#' * pheno: (data.frame) phenotype matrix as used by predictor
#' * TT_STATUS: initial train/test split assignment
#' * netScores: (list) keys are classes, value is a data.frame with 
#' networks and scores following feature selection. Nets with score>0 are
#' shown.
#' * confmat: (matrix) confusion matrix with one row per cutoff
#' * predRanks: (list) GM rankings of test sampels for each network cutoff
#' In addition, all the intervening work will be stored in <outDir>
#' @export
buildPredictor_resampling <- function(pheno, pdat, predClass, unitSets,
	pctT=0.7,numResamples=3L, nFoldCV=10L,numCores=1L,GMmemory=4L,
	outDir=".",overwrite=FALSE,...) {
DEBUG_FLAG <- TRUE

if (file.exists(outDir)) {
	if (!overwrite) stop("output directory exists. Choices: provide path to non-existing directory, set overwrite option to TRUE, or manually delete this directory")
	else unlink(outDir,recursive=TRUE)
}
dir.create(outDir)

pheno$STATUS[which(!pheno$STATUS %in% predClass)] <- "other"
subtypes <- c(predClass,"other")

pheno$TT_STATUS <- splitTestTrain(pheno,
    	pctT = pctT,setSeed=42,predClass=predClass)

pheno_FULL	<- pheno
pdat_FULL 	<- pdat
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
pdat		<- pdat[,which(colnames(pdat)%in% pheno$ID)]

# ############################
# Feature selection

# run feature selection with resamplings, once per class
featNets <- list()
resampTT <- list()
for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g)
    	if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
	dir.create(pDir)

	cat(sprintf("\n ********** \nClass: %s\n ********** \n",g))
	TT_STATUS <- splitTestTrain_partition(
		pheno,nFold=numResamples, predClass=g,setSeed=103)
	save(TT_STATUS, file=sprintf("%s/TT_STATUS_resampling.Rdata",pDir))
	resampTT[[g]] <- TT_STATUS
	
	for (repNum in 1:numResamples) {
		cat(sprintf("\nResampling %i ***********\n",repNum))
		pheno_subtype <- pheno
		pDir <- sprintf("%s/%s/part%i", outDir,g,repNum)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
		dir.create(pDir)
	
		pheno_subtype$INNER_TT_STATUS <- TT_STATUS[[repNum]]
		pheno_subtype <- subset(pheno_subtype, 
					INNER_TT_STATUS%in% "TRAIN")
		cat(sprintf("\t%i training\n",nrow(pheno_subtype)))
		
		## label patients not in the current class as a residual
		tmp <- which(!pheno_subtype$STATUS %in% g)
		pheno_subtype$STATUS[tmp]<-"nonpred"
		print(table(pheno_subtype$STATUS,useNA="always"))
	
		# prepare nets
		profDir <- sprintf("%s/profiles",pDir)
		netDir <- sprintf("%s/networks",pDir)
		pdat_tmp	<- pdat[,which(colnames(pdat)%in% pheno_subtype$ID)]
		netList <- makePSN_NamedMatrix(pdat_tmp, rownames(pdat_tmp), 
		        unitSets,profDir,verbose=FALSE,
		        numCores=numCores,writeProfiles=TRUE,...)
		netList <- unlist(netList)
	
		dbDir	<- GM_createDB(profDir, pheno_subtype$ID, outDir,
							numCores=numCores)
		resDir	<- sprintf("%s/GM_results",pDir)

		## cross validation using training samples from predClass
		isg <- which(pheno_subtype$STATUS %in% g)
		trainPred <- pheno_subtype$ID[isg]
		GM_runCV_featureSet(trainPred, resDir, dbDir$dbDir, 
			nrow(pheno_subtype),verbose=T, numCores=numCores,
			GMmemory=GMmemory,nFold=nFoldCV)
		
	    # Compute network score
		nrank <- dir(path=resDir,pattern="NRANK$")
		pTally	<- GM_networkTally(paste(resDir,nrank,sep="/"))
		# write to file
		tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
		write.table(pTally,file=tallyFile,sep="\t",col=T,row=F,quote=F)
	}
	
	## tally total network score across resamplings
	pathScore <- list()
	for (currRep in 1:numResamples) {
		f <- sprintf("%s/%s/part%i/GM_results/%s_pathway_CV_score.txt",
			outDir,g,currRep,g)
		dat <- read.delim(f,sep="\t",header=TRUE,as.is=TRUE)
		for (k in 1:nrow(dat)) {
			p <- dat[k,1]
			if (p %in% names(pathScore)) 
				pathScore[[p]] <-pathScore[[p]]+dat[k,2]
			else pathScore[[p]] <- dat[k,2]
		}
	}
	df <- data.frame(PATHWAY_NAME=names(pathScore), SCORE=unlist(pathScore))
	featNets[[g]] <- df
	write.table(df,file=sprintf("%s/%s/%s_pathwayScore.txt", outDir,
			g,g),sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}
cat("* Feature selection complete\n")

# ############################
# Class prediction

# now create GM databases for each class
# should contain train + test patients
# and be limited to nets that pass feature selection at all thresholds
pheno <- pheno_FULL
predRes <- list() 	## predRes[[cutoff]] contains predictions for 
					## score >= cutoff.
maxScore <- numResamples * nFoldCV
for (cutoff in 1:maxScore) predRes[[cutoff]] <- list()
# get patient rankings at each cutoff
for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g)
	pTally <- read.delim(
		sprintf("%s/%s_pathwayScore.txt",pDir,g),
		sep="\t",h=T,as.is=T)
	pTally <- pTally[which(pTally[,2]>=1),]
	cat(sprintf("%s: %i pathways\n",g,nrow(pTally)))

	# create new db
	profDir <- sprintf("%s/profiles",pDir)
	cleanName <- sub(".profile","",pTally[,1])
	cleanName <- sub("_cont","",cleanName)
	tmp <- makePSN_NamedMatrix(pdat_FULL,rownames(pdat),
		unitSets[which(names(unitSets)%in% cleanName)],
		profDir,verbose=F,numCores=numCores,writeProfiles=TRUE,...)
	dbDir <- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

	# query of all training samples for this class
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
							 pheno$TT_STATUS%in%"TRAIN")]
	
	cl <- makeCluster(numCores)
	registerDoParallel(cl)
	foreach(cutoff=1:maxScore, .packages="netDx") %dopar% {
		cat(sprintf("\tCutoff = %i\n",cutoff))
		qFile <- sprintf("%s/%s_cutoff%i",pDir,g,cutoff)
		curr_p <- pTally[which(pTally[,2]>=cutoff),1]
		if (length(curr_p)>0){
			netDx::GM_writeQueryFile(qSamps,curr_p,nrow(pheno),qFile)
			resFile <- netDx::runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
			system(sprintf("unlink %s", resFile))
		}
	}
	stopCluster(cl)

	# collect rankings
	for (cutoff in 1:maxScore) {
		resFile <- sprintf("%s/%s_cutoff%i-results.report.txt.PRANK",
			pDir,g,cutoff)
		predRes[[cutoff]][[g]] <- GM_getQueryROC(resFile,pheno,g)
	}
}

save(predRes, pheno,file=sprintf("%s/predictionResults.Rdata",outDir))

# compile performance at each cutoff
outClass <- list()
outmat <- matrix(NA, nrow=maxScore, ncol=9)
colnames(outmat) <- c("score","tp","tn","fp","fn","tpr","fpr",
					  "accuracy","ppv")
for (cutoff in 1:maxScore) {
	outClass[[cutoff]] <- GM_OneVAll_getClass(predRes[[cutoff]])
	both <- merge(x=pheno,y=outClass[[cutoff]],by="ID")
	print(table(both[,c("STATUS","PRED_CLASS")]))
	pos <- (both$STATUS %in% predClass)
	tp <- sum(both$PRED_CLASS[pos]==predClass)
	fp <- sum(both$PRED_CLASS[!pos]==predClass)
	tn <- sum(both$PRED_CLASS[!pos]=="other")
	fn <- sum(both$PRED_CLASS[pos]=="other")
	acc <- (tp+tn)/nrow(both)
	ppv <- tp/(tp+fp)
	outmat[cutoff,] <- c(cutoff,tp,tn,fp,fn,tp/(tp+fn), fp/(fp+tn),acc,ppv)
}

cat("* Predictions complete!\n")
out <- list(pheno=pheno_FULL,TT_STATUS=resampTT,netScore=featNets,
		confmat=outmat,predRanks=outClass)

out
}

