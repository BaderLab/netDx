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
#' 		compile features into database for cross-validation
#'		score networks out of 10
#' end
#' using test samples from all resamplings, measure predictor performance.
#' 
#' In short, this function performs all steps involved in building and 
#' evaluating the predictor. 
#' @param phenoDF (data.frame) sample metadat. patient ID,STATUS
#' @param cnv_GR (GRanges) genetic events. Must contain "ID" column mapping
#' the event to a patient. ID must correspond to the ID column in phenoDF
#' @param predClass (char) patient class to predict
#' @param outDir (char) path to dir where results should be stored. 
#' Results for resampling i are under \code{<outDir>/part<i>}, while
#' predictor evaluation results are directly in \code{outDir}.
#' @param splitN (integer) number of data resamplings to use
#' @param featScoreMax (integer) max score for features in feature selection
#' @param filter_WtSum (numeric between 5-100) Limit to top-ranked 
#' networks such that cumulative weight is less than this parameter. 
#' e.g. If filter_WtSum=20, first order networks by decreasing weight; 
#' then keep those whose cumulative weight <= 20.
#' @param enrichLabels (logical) if TRUE, applies label enrichment to train
#' networks
#' @param enrichPthresh (numeric between 0 and 1) networks with label
#' enrichment p-value below this threshold  pass enrichment
#' @param numPermsEnrich (integer) number of permutations for label enrichment
#' @param minEnr (integer -1 to 1) minEnr param in enrichLabelsNets()
#' @param numCores (integer) num cores for parallel processing
#' @param FS_numCores (integer) num cores for running GM. If NULL, is set
#' to max(1,numCores-1). Set to a lower value if the default setting
#' gives out-of-memory error. This may happen if networks are denser than
#' expected
#' @param ... params for runFeatureSelection()
#' @return No value. Side efect of writing feature scores and predictions to 
#' `outDir`. Output format is that of running `RR_featureTally()`. Refer to 
#' that function for the details of output files.
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @export
buildPredictor_sparseGenetic <- function(phenoDF,cnv_GR,predClass,
	group_GRList,outDir=tempdir(),
	splitN=3L, featScoreMax=10L,
	filter_WtSum=100L,
	enrichLabels=TRUE,enrichPthresh=0.07,numPermsEnrich=2500L,minEnr=-1,
	numCores=1L,FS_numCores=NULL,...) {

	netDir <- sprintf("%s/networks_orig",outDir)
	netList <- makePSN_RangeSets(cnv_GR, group_GRList,netDir,verbose=FALSE)

	p 	<- countPatientsInNet(netDir,netList, phenoDF$ID)
	tmp	<- updateNets(p,phenoDF,writeNewNets=FALSE,verbose=FALSE)

	netmat	<- tmp[[1]]
	phenoDF	<- tmp[[2]] 

	if (is.null(FS_numCores)) FS_numCores <- max(1,numCores-1)
	
	# split into testing and training - resampling mode
	message("* Resampling train/test samples")
	TT_STATUS 	<- splitTestTrain_resampling(phenoDF, nFold=splitN,
		predClass=predClass, verbose=TRUE)
	p_full 		<- netmat
	pheno_full	<- phenoDF
	if (any(colnames(pheno_full)%in% "TT_STATUS")) {
		message(paste("** Warning, found TT_STATUS column. ",
			"netDx adds its own column so this one will be removed **",sep=""))
		pheno_full <- pheno_full[,-which(colnames(pheno_full)%in%"TT_STATUS")]
	}

	pScore		<- list()
	enrichedNets	<- list()
	for (k in seq_len(length(TT_STATUS))) {
		p 		<- p_full
		pheno 	<- pheno_full
	
		pheno  <- cbind(pheno, TT_STATUS=TT_STATUS[[k]])
		message("----------------------------------------")
		message(sprintf("Resampling round %i", k))
		message("----------------------------------------")
		print(table(pheno[,c("STATUS","TT_STATUS")]))
	
		newOut <- sprintf("%s/part%i",outDir,k)
		dir.create(newOut)

		# write patient status for this round. 
		outF <- sprintf("%s/TT_STATUS.txt",newOut)
		write.table(pheno,file=outF,sep="\t",col.names=TRUE,
			row.names=FALSE,quote=FALSE)
	
		message("# patients: train only")
		pheno_train <- subset(pheno, TT_STATUS %in% "TRAIN")
		p_train  	<- p[which(rownames(p) %in% pheno_train$ID),]
		print(nrow(pheno_train))
		
		# update nets
		message("Training only:")
		trainNetDir <- sprintf("%s/networks",newOut)
		tmp			<- updateNets(p_train,pheno_train, 
							oldNetDir=netDir, newNetDir=trainNetDir,verbose=FALSE)
		p_train		<- tmp[[1]]
		pheno_train	<- tmp[[2]]

		# label enrichment
		if (enrichLabels) {
			message("Running label enrichment")
			tmpDir <- sprintf("%s/tmp",outDir)
			if (!file.exists(tmpDir)) dir.create(tmpDir)
			netInfo <- enrichLabelNets(trainNetDir,pheno_train,newOut,
				predClass=predClass,numReps=numPermsEnrich,
				numCores=numCores,tmpDir=tmpDir,verbose=FALSE)
			pvals   <- as.numeric(netInfo[,"pctl"])

			netInfo <- netInfo[which(pvals < enrichPthresh),] 
			print(nrow(netInfo))
			p_train  	<- p_train[,
				which(colnames(p_train) %in% rownames(netInfo))]
		
			# update nets after enrichment
			trainNetDir <- sprintf("%s/networksEnriched",newOut)
			tmp			<- updateNets(p_train, pheno_train,
							oldNetDir=netDir, 
							newNetDir=trainNetDir,verbose=FALSE)
			p_train		<- tmp[[1]]
			pheno_train	<- tmp[[2]]

			enrichedNets[[k]]	<- rownames(netInfo)
		}

		pheno_train <- setupFeatureDB(pheno_train,trainNetDir)
		moveInteractionNets(netDir=trainNetDir,
				outDir=sprintf("%s/INTERACTIONS",trainNetDir),
				pheno=pheno_train)
		
		# create networks for cross-validation
		x <- compileFeatures(trainNetDir,newOut,verbose=FALSE)
		
		# we query for training samples of the predictor class
		trainPred <- pheno_train$ID[
			which(pheno_train$STATUS %in% predClass)]
		resDir    <- sprintf("%s/GM_results",newOut)
		dbPath     <- sprintf("%s/dataset",newOut)
		t0 <- Sys.time()
		runFeatureSelection(trainID_pred=trainPred, 
				outDir=resDir, dbPath=dbPath, 
				numTrainSamps=nrow(p_train),
				verbose=FALSE,numCores=FS_numCores,
				featScoreMax=featScoreMax,...)
		t1 <- Sys.time()
		message("Score features for this train/test split")
		print(t1-t0)
		
		# collect results
		nrankFiles	<- paste(resDir,dir(path=resDir,pattern="NRANK$"),
			sep="/")
		pathwayRank	<- compileFeatureScores(nrankFiles,
			filter_WtSum=filter_WtSum,verbose=FALSE)
		write.table(pathwayRank,file=sprintf("%s/pathwayScore.txt",resDir),
			col.names=TRUE,row.names=FALSE,quote=FALSE)
		
		pScore[[k]]	<- pathwayRank
	}

	phenoDF$TT_STATUS <- NULL
	#outDat <- sprintf("%s/resampling_savedDat.Rdata",outDir)

	out <- list(
		netmat=netmat,
		pheno=phenoDF,
		TT_STATUS=TT_STATUS,
		pathwayScores=pScore,
		enrichedNets=enrichedNets)	
	#save(netmat,phenoDF,TT_STATUS,predClass,pScore,outDir,enrichLabels,
	#	 enrichedNets, file=outDat)
# Measure performance based on pathway tally across pplits
	RR_featureTally(netmat, phenoDF, TT_STATUS, predClass,pScore,
			outDir,enrichLabels, enrichedNets,verbose=FALSE)
}
