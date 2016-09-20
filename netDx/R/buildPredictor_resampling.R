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
#' @param seed_trainTest (integer) RNG seed for initial train/test split.
#' Set to NULL to use current state of RNG. Set to a constant for 
#' reproducibility
#' @param seed_resampling (integer) RNG seed for training assignment for
#' of the resamplings (splitTrainTest_partition() method)
#' @param seed_CVqueries (integer) RNG seed for queries in k-fold cross
#' validation (used by makeCVqueries() which is called by 
#' GM_runCV_featureSet())
#' @param ... parameters for makePSN_NamedMatrix()
#' @return (list) Includes:
#' * pheno: (data.frame) phenotype matrix as used by predictor
#' * TT_STATUS: initial train/test split assignment
#' * netScores: (list) keys are classes, value is a data.frame with 
#' networks and scores following feature selection. Nets with score>0 are
#' shown.
#' * confmat: (matrix) confusion matrix with one row per cutoff
#' * predRanks: (list) GM rankings of test samples for each network cutoff
#' In addition, all the intervening work will be stored in <outDir>
#' @export
buildPredictor_resampling <- function(pheno, pdat, predClass, unitSets,
	pctT=0.7,numResamples=3L, nFoldCV=10L,numCores=1L,GMmemory=4L,
	outDir=".",overwrite=FALSE,seed_trainTest=42L,seed_resampling=103L,
	seed_CVqueries=42L,...) {
DEBUG_FLAG <- TRUE

if (file.exists(outDir)) {
	if (!overwrite) stop("output directory exists. Choices: provide path to non-existing directory, set overwrite option to TRUE, or manually delete this directory")
	else unlink(outDir,recursive=TRUE)
}
dir.create(outDir)

pheno$STATUS[which(!pheno$STATUS %in% predClass)] <- "other"
subtypes <- c(predClass,"other")

pheno$TT_STATUS <- splitTestTrain(pheno,
    	pctT = pctT,setSeed=seed_trainTest,predClass=predClass)

pheno_FULL	<- pheno
pdat_FULL 	<- pdat
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
pdat		<- pdat[,which(colnames(pdat)%in% pheno$ID)]

# --------------------------------------------------
# Phase 1. Feature selection, assigning net scores
# --------------------------------------------------

# run feature selection with resamplings, once per class
featNets <- list()
resampTT <- splitTestTrain_partition(
	pheno,nFold=numResamples, predClass=predClass,setSeed=seed_resampling)
save(resampTT, file=sprintf("%s/TT_STATUS_resampling.Rdata",outDir))
for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g)
    	if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
	dir.create(pDir)

	cat(sprintf("\n ********** \nClass: %s\n ********** \n",g))
	for (repNum in 1:numResamples) {
		cat(sprintf("\nResampling %i ***********\n",repNum))
		pheno_subtype <- pheno
		pDir <- sprintf("%s/%s/part%i", outDir,g,repNum)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
		dir.create(pDir)
	
		pheno_subtype$INNER_TT_STATUS <- resampTT[[repNum]]
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
			GMmemory=GMmemory,nFold=nFoldCV,setSeed=seed_CVqueries)
		
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
	write.table(df,file=sprintf("%s/%s_pathwayScore.txt", outDir,
			g),sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}

##  pathways by subtype
netScores <- list()
for (g in subtypes) {
	pTally <- read.delim(
		sprintf("%s/%s_pathwayScore.txt",outDir,g),
		sep="\t",h=T,as.is=T)
	pTally <- pTally[which(pTally[,2]>=1),]
	pTally[,1] <- sub(".profile","",pTally[,1])
	netScores[[g]] <- pTally
}

cat("* Net scores computed\n")

# --------------------------------------------------
# Phase 2. Cutoff selection from average of resampling "test"
# --------------------------------------------------
maxScore <- numResamples * nFoldCV
for (g in names(netScores)) 
		maxScore <- min(maxScore,max(netScores[[g]][,2]))
cat(sprintf("maxScore = %i\n", maxScore))
perf_resamp <- list()
for (rep in 1:numResamples) {
	perfDir <- sprintf("%s/eval/part%i",outDir,rep)
	if (file.exists(perfDir)) unlink(perfDir,recursive=TRUE)
	dir.create(perfDir,recursive=TRUE)
	pheno_cur <- pheno
	pheno_cur$TT_STATUS <- resampTT[[rep]]
	print(table(pheno_cur[,c("STATUS","TT_STATUS")]))
	perf_resamp[[rep]] <- GM_predClass_cutoffs(pheno_cur,pdat,
		predClass=predClass,netScores=netScores,unitSets=unitSets,
		maxScore=maxScore,outDir=perfDir,numCores=numCores)
}
save(perf_resamp,file=sprintf("%s/resamplingPerf.Rdata",outDir))

# pick cutoff with best mean performance across resamplings
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

# --------------------------------------------------
# Phase 3. Model evaluation on the blind test samples
# --------------------------------------------------
cat("\n\n* Model evaluation\n")
finalDir <- sprintf("%s/test", outDir)
finalNets <- lapply(netScores, function(x) x[which(x[,2]>=bestCutoff),1])
testRes <- GM_predClass_once(pheno_FULL,pdat_FULL,predClass,unitSets,
	patNets=finalNets,outDir=finalDir,numCores=numCores)

cat("* Predictions complete!\n")
out <- list(pheno=pheno_FULL,resampTT=resampTT,netScore=featNets,
		perfResamp=perf_resamp, bestCutoff=bestCutoff, testRes=testRes)
		
out
}

