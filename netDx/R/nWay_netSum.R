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
#' @param seed_resampling (integer) RNG seed for deciding hold-out sets
#' while resampling.
#' @param nFoldCV (integer) number of folds in the inner cross-validation
#' loop
#' @param filter_WtSum (numeric between 5-100) Limit to top-ranked 
#' networks such that cumulative weight is less than this parameter. 
#' e.g. If filter_WtSum=20, first order networks by decreasing weight; 
#' then keep those whose cumulative weight <= 20.
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
#' @param seed_CVqueries (integer) RNG seed for inner cross-validation
#' loop
#' @param ... params for runFeatureSelection()
#' @importFrom reshape2 melt
#' @export
Nway_netSum <- function(netmat=NULL, phenoDF,predClass,outDir,netDir,
	splitN=3L,seed_resampling=103L, nFoldCV=10L,filter_WtSum=100L,
	cliqueFilter=TRUE,cliquePthresh=0.07,cliqueReps=2500L,minEnr=-1,
	numCores=1L,GM_numCores=NULL,useAttributes=NULL,
	seed_CVqueries=42L,...) {

	if (is.null(GM_numCores)) GM_numCores <- max(1,numCores-1)
	
	# split into testing and training - resampling mode
	cat("* Resampling train/test samples\n")
	TT_STATUS 	<- splitTestTrain_resampling(phenoDF, nFold=splitN,
		predClass=predClass, setSeed=seed_resampling, verbose=TRUE)
	p_full 		<- netmat
	pheno_full	<- phenoDF
	if (any(colnames(pheno_full)%in% "TT_STATUS")) {
		cat("** Warning, found TT_STATUS column. netDx adds its own column so this one will be removed **\n")
		pheno_full <- pheno_full[,-which(colnames(pheno_full)%in%"TT_STATUS")]
	}

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
			tmpDir <- sprintf("%s/tmp",outDir)
			if (!file.exists(tmpDir)) dir.create(tmpDir)
			netInfo <- cliqueFilterNets(trainNetDir,pheno_train,newOut,
				predClass=predClass,numReps=cliqueReps,numCores=numCores,
				tmpDir=tmpDir)
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
		x <- compileFeatures(trainNetDir, rownames(p_train), newOut)
					
		
		# we query for training samples of the predictor class
		trainPred <- pheno_train$ID[
			which(pheno_train$STATUS %in% predClass)]
		resDir    <- sprintf("%s/GM_results",newOut)
		dbPath     <- sprintf("%s/dataset",newOut)
		t0 <- Sys.time()
		runFeatureSelection(trainPred, resDir, dbPath, 
				nrow(p_train),verbose=TRUE,numCores=GM_numCores,
				nFold=nFoldCV,seed_CVqueries=seed_CVqueries,...)
		t1 <- Sys.time()
		cat("Time to run inner CV loop:\n")
		print(t1-t0)
		
		# collect results
		nrankFiles	<- paste(resDir,dir(path=resDir,pattern="NRANK$"),
			sep="/")
		pathwayRank	<- compileFeatureScores(nrankFiles,
			filter_WtSum=filter_WtSum,verbose=TRUE)
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
			outDir,cliqueFilter, cliqueNets,verbose=FALSE)
}
