#' Performs feature selection using multiple resamplings of the data
#'
#' @details This function is used for feature selection of patient networks,
#' using multiple resamplings of input data. It is intended for use in
#' the scenario where patient networks are sparse and binary. 
#' This function should be called after defining all patient networks. It
#' performs the following steps:
#' For i = 1..splitN
#' 		randomly split patients into training and test
#' 		(optional) filter training networks to exclude random-like networks
#' 		create GeneMANIA database for cross-validation
#' 		run 10-fold cross validation
#'		score networks based on 10-fold CV
#' end
#' using test samples from all resamplings, measure predictor performance.
#' 
#' In short, this function performs all steps involved in building and 
#' evaluating the predictor. 
#' @param netmat (matrix) output of countPatientsInNet()
#' @param phenoDF (data.frame) patient ID,STATUS, ATTRIB_GROUP,ATTRIB_NAME
#' last two are optional and only used if useAttributes is not NULL. Note
#' that duplicate ATTRIB_NAMEs are currently not supported, even if they
#' are in different ATTRIB_GROUPs.
#' @param predClass (char) patient class to predict
#' @param outDir (char) path to dir where results should be stored. 
#' Results for resampling i are under \code{<outDir>/part<i>}, while
#' predictor evaluation results are directly in \code{outDir}.
#' @param netDir (char) path to dir where patient similarity networks
#' are stored. These networks should contain all training and test patients
#' in the dataset. This function will subset these networks to generate
#' networks consisting e.g. only of training samples.
#' @param splitN (integer) number of data resamplings to use
#' @param nFoldCV (integer) number of folds in the inner cross-validation
#' loop
#' @param cliqueFilter (logical) if TRUE, applies clique filtering to train
#' networks
#' @param cliquePthresh (numeric between 0 and 1) networks with clique-
#' filtering p-value below this threshold  pass clique-filtering
#' @param cliqueReps (integer) number of permutations for clique filtering
#' @param minEnr (integer -1 to 1) minEnr param in cliqueFilterNets()
#' @param numCores (integer) num cores for parallel processing
#' @param GM_numCores (integer) num cores for running GM. If NULL, is set
#' to max(1,numCores-1). Set to a lower value if the default setting
#' gives out-of-memory error. This may happen if networks are denser than
#' expected
#' @param useAttributes (char) vector of attribute names to be used in 
#' executing GM queries. Note: Not currently well-tested, suggest leaving
#' as NULL.
#' @param ... params for GM_runCV_featureSet()
#' @importFrom reshape2 melt
#' @export
Nway_netSum <- function(netmat=NULL, phenoDF,predClass,outDir,netDir,
	splitN=3L,nFoldCV=10L,
	cliqueFilter=TRUE,cliquePthresh=0.07,cliqueReps=2500L,minEnr=-1,
	numCores=1L,GM_numCores=NULL,useAttributes=NULL,...) {

		if (is.null(GM_numCores)) GM_numCores <- max(1,numCores-1)
	
	# split into testing and training - resampling mode
	cat("* Resampling train/test samples\n")
	TT_STATUS 	<- splitTestTrain_partition(phenoDF, nFold=splitN,
		predClass=predClass, verbose=TRUE)
	p_full 		<- netmat
	pheno_full	<- phenoDF
	pheno_full <- pheno_full[,-which(colnames(pheno_full)%in%"TT_STATUS")]

	pScore		<- list()
	cliqueNets	<- list()
	for (k in 1:length(TT_STATUS)) {
		p 		<- p_full
		pheno 	<- pheno_full
	
		pheno  <- cbind(pheno, TT_STATUS=TT_STATUS[[k]])
		cat("----------------------------------------\n")
		cat(sprintf("Resampling round %i\n", k))
		cat("----------------------------------------\n")
		print(table(pheno[,c("STATUS","TT_STATUS")]))
	
		newOut <- sprintf("%s/part%i",outDir,k)
		dir.create(newOut)

		# write patient status for this round. 
		outF <- sprintf("%s/TT_STATUS.txt",newOut)
		write.table(pheno,file=outF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
	
		cat("# patients: train only\n")
		pheno_train <- subset(pheno, TT_STATUS %in% "TRAIN")
		p_train  	<- p[which(rownames(p) %in% pheno_train$ID),]
		print(nrow(pheno_train))
		
		# update nets
		cat("Training only:\n")
		trainNetDir <- sprintf("%s/networks",newOut)
		tmp			<- updateNets(p_train,pheno_train, 
							oldNetDir=netDir, newNetDir=trainNetDir)
		p_train		<- tmp[[1]]
		pheno_train	<- tmp[[2]]
		
		## ----------------------------------------------------------------
		# clique filtering
		if (cliqueFilter) {
			cat("Running clique-filtering\n")
			netInfo <- cliqueFilterNets(trainNetDir,pheno_train,newOut,
				predClass=predClass,numReps=cliqueReps,numCores=numCores)
			pvals   <- as.numeric(netInfo[,"pctl"])

			netInfo <- netInfo[which(pvals < cliquePthresh),] 
			print(nrow(netInfo))
			p_train  	<- p_train[,
				which(colnames(p_train) %in% rownames(netInfo))]
		
			# update nets after clique-filtering
			cat("Clique filtered\n")
			trainNetDir <- sprintf("%s/networksCliqueFilt",newOut)
			tmp			<- updateNets(p_train, pheno_train,
							oldNetDir=netDir, newNetDir=trainNetDir)
			p_train		<- tmp[[1]]
			pheno_train	<- tmp[[2]]

			cliqueNets[[k]]	<- rownames(netInfo)
		}

		# create networks for cross-validation
		attrib_DF <- NULL
		if (!is.null(useAttributes)) {
				cat("not yet implemented")
				browser()
		}
		x <- GM_createDB(trainNetDir, rownames(p_train), newOut)
					
		
		# we query for training samples of the predictor class
		trainPred <- pheno_train$ID[
			which(pheno_train$STATUS %in% predClass)]
		resDir    <- sprintf("%s/GM_results",newOut)
		GM_db     <- sprintf("%s/dataset",newOut)
		t0 <- Sys.time()
		GM_runCV_featureSet(trainPred, resDir, GM_db, 
				nrow(p_train),verbose=TRUE,numCores=GM_numCores,
				nFold=nFoldCV,...)
		t1 <- Sys.time()
		cat("Time to run inner CV loop:\n")
		print(t1-t0)
		
		# collect results
		nrankFiles	<- paste(resDir,dir(path=resDir,pattern="NRANK$"),
			sep="/")
		pathwayRank	<- GM_networkTally(nrankFiles)
		write.table(pathwayRank,file=sprintf("%s/pathwayScore.txt",resDir),
			col=T,row=F,quote=F)
		
		pScore[[k]]	<- pathwayRank
	}

	phenoDF$TT_STATUS <- NULL
	outDat <- sprintf("%s/resampling_savedDat.Rdata",outDir)
	save(netmat,phenoDF,TT_STATUS,predClass,pScore,outDir,cliqueFilter,
		 cliqueNets, file=outDat)
# Measure performance based on pathway tally across splits
	RR_featureTally(netmat, phenoDF, TT_STATUS, predClass,pScore,
			outDir,cliqueFilter, cliqueNets)
}
