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
#' @param featScoreMax (integer) max score for features in feature selection
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
#' @param FS_numCores (integer) num cores for running GM. If NULL, is set
#' to max(1,numCores-1). Set to a lower value if the default setting
#' gives out-of-memory error. This may happen if networks are denser than
#' expected
#' @param seed_queryResample (integer) RNG seed for inner cross-validation
#' loop
#' @param ... params for runFeatureSelection()
#' @return No value. Side efect of writing feature scores and predictions to 
#' `outDir`. Output format is that of running `RR_featureTally()`. Refer to 
#' that function for the details of output files.
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @export
Nway_netSum <- function(netmat=NULL, phenoDF,predClass,
	outDir=tempdir(),netDir,
	splitN=3L,seed_resampling=103L, featScoreMax=10L,
	filter_WtSum=100L,
	cliqueFilter=TRUE,cliquePthresh=0.07,cliqueReps=2500L,minEnr=-1,
	numCores=1L,FS_numCores=NULL,
	seed_queryResample=42L,...) {

	if (is.null(FS_numCores)) FS_numCores <- max(1,numCores-1)
	
	# split into testing and training - resampling mode
	message("* Resampling train/test samples")
	TT_STATUS 	<- splitTestTrain_resampling(phenoDF, nFold=splitN,
		predClass=predClass, verbose=TRUE)
	p_full 		<- netmat
	pheno_full	<- phenoDF
	if (any(colnames(pheno_full)%in% "TT_STATUS")) {
		message("** Warning, found TT_STATUS column. netDx adds its own column so this one will be removed **")
		pheno_full <- pheno_full[,-which(colnames(pheno_full)%in%"TT_STATUS")]
	}

	pScore		<- list()
	cliqueNets	<- list()
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
							oldNetDir=netDir, newNetDir=trainNetDir)
		p_train		<- tmp[[1]]
		pheno_train	<- tmp[[2]]
		

		# clique filtering
		if (cliqueFilter) {
			message("Running clique-filtering")
			tmpDir <- sprintf("%s/tmp",outDir)
			if (!file.exists(tmpDir)) dir.create(tmpDir)
			netInfo <- cliqueFilterNets(trainNetDir,pheno_train,newOut,
				predClass=predClass,numReps=cliqueReps,
				numCores=numCores,tmpDir=tmpDir)
			pvals   <- as.numeric(netInfo[,"pctl"])

			netInfo <- netInfo[which(pvals < cliquePthresh),] 
			print(nrow(netInfo))
			p_train  	<- p_train[,
				which(colnames(p_train) %in% rownames(netInfo))]
		
			# update nets after clique-filtering
			trainNetDir <- sprintf("%s/networksCliqueFilt",newOut)
			tmp			<- updateNets(p_train, pheno_train,
							oldNetDir=netDir, 
							newNetDir=trainNetDir)
			p_train		<- tmp[[1]]
			pheno_train	<- tmp[[2]]

			cliqueNets[[k]]	<- rownames(netInfo)
		}

		pheno_train <- setupFeatureDB(pheno_train,trainNetDir)
		moveInteractionNets(netDir=trainNetDir,
				outDir=sprintf("%s/INTERACTIONS",trainNetDir),
				pheno=pheno_train)
		
		# create networks for cross-validation
		x <- compileFeatures(trainNetDir,newOut)
message("compile done")
		
		# we query for training samples of the predictor class
		trainPred <- pheno_train$ID[
			which(pheno_train$STATUS %in% predClass)]
		resDir    <- sprintf("%s/GM_results",newOut)
		dbPath     <- sprintf("%s/dataset",newOut)
		t0 <- Sys.time()
		runFeatureSelection(trainID_pred=trainPred, 
				outDir=resDir, dbPath=dbPath, 
				numTrainSamps=nrow(p_train),
				verbose=TRUE,numCores=FS_numCores,
				featScoreMax=featScoreMax,...)
		t1 <- Sys.time()
		message("Score features for this train/test split")
		print(t1-t0)
		
		# collect results
		nrankFiles	<- paste(resDir,dir(path=resDir,pattern="NRANK$"),
			sep="/")
		pathwayRank	<- compileFeatureScores(nrankFiles,
			filter_WtSum=filter_WtSum,verbose=TRUE)
		write.table(pathwayRank,file=sprintf("%s/pathwayScore.txt",resDir),
			col.names=TRUE,row.names=FALSE,quote=FALSE)
		
		pScore[[k]]	<- pathwayRank
	}

	phenoDF$TT_STATUS <- NULL
	outDat <- sprintf("%s/resampling_savedDat.Rdata",outDir)
	save(netmat,phenoDF,TT_STATUS,predClass,pScore,outDir,cliqueFilter,
		 cliqueNets, file=outDat)
# Measure performance based on pathway tally across pplits
	RR_featureTally(netmat, phenoDF, TT_STATUS, predClass,pScore,
			outDir,cliqueFilter, cliqueNets,verbose=FALSE)
}
