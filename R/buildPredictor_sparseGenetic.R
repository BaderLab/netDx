#' Performs feature selection using multiple resamplings of the data
#'
#' @details This function is used for feature selection of patient networks,
#' using multiple resamplings of input data. It is intended for use in
#' the scenario where patient networks are sparse and binary. 
#' This function should be called after defining all patient networks. It
#' performs the following steps:
#' For i = 1..numSplits
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
#' @param group_GRList (list) List of GRangesList indicating grouping 
#' rules for CNVs.
#' For example, in a pathway-based design, each key value would be a pathway
#' name, and the value would be a RangesList containing coordinates of the
#' member genes
#' @param predClass (char) patient class to predict
#' @param outDir (char) path to dir where results should be stored. 
#' Results for resampling i are under \code{<outDir>/part<i>}, while
#' predictor evaluation results are directly in \code{outDir}.
#' @param numSplits (integer) number of data resamplings to use
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
#' @return (list) Predictor results
#' 1) phenoDF (data.frame): subset of phenoDF provided as input, but limited
#' to patients that have at least one event in the universe of possibilities
#' e.g. if using pathway-level features, then this table excludes patients
#' with zero CNVs in pathways
#' 2) netmat (data.frame): Count of genetic events by patients (rows) in
#' pathways (columns). Used as input to the feature selection algorithm
#' 3) pathwayScores (list): pathway scores for each of the data splits. 
#' Each value in the list is a data.frame containing pathway names and
#' scores.
#' 4) enrichedNets (list): This entry is only found if enrichLabels is set
#' to TRUE. It contains the vector of features that passed label enrichment
#' in each split of the data. 
#' 5 - 9) Output of RR_featureTally:
#' 5) cumulativeFeatScores: pathway name, cumulative score over 
#' N-way data resampling.
#' 6) performance_denAllNets: positive,negative calls at each cutoff: 
#' network score cutoff (score); num networks at cutoff (numPathways) ; 
#' total +, ground truth (pred_tot); + calls (pred_ol); 
#' + calls as pct of total (pred_pct); total -, ground truth (other_tot) ; 
#' - calls  (other_ol) ; - calls as pct of total (other_pct) ; ratio of 
#' pred_pct and other_pct (rr) ; min. pred_pct in all resamplings 
#' (pred_pct_min) ; max pred_pct in all resamplings (pred_pct_max) ; 
#' min other_pct in all resamplings (other_pct_min); max other_pct in all
#' resamplings (other_pct_max)
#' 7) performance_denEnrichedNets: positive, negative calls at each cutoff 
#' label enrichment option: format same as performance_denAllNets. 
#' However, the denominator here is limited to
#' patients present in networks that pass label enrichment
#' 8) resamplingPerformance: breakdown of performance for each of the 
#' resamplings, at each of the cutoffs. 
#' This is a list of length 2, one for allNets and one for enrichedNets. 
#' The value is a matrix with (resamp * 7) columns and S rows,
#' one row per score. The columns contain the following information
#' per resampling:
#'	 1) pred_total: total num patients of predClass
#'	 2) pred_OL: num of pred_total with a CNV in the selected net
#'	 3) pred_OL_pct: 2) divided by 1) (percent)
#'	 4) other_total: total num patients of other class(non-predClass)
#'	 5) other_OL: num of other_total with CNV in selected net
#'	 6) other_OL_pct: 5) divided by 4) (percent)
#'	 7) relEnr: 6) divided by 3).
#' 
#' @importFrom reshape2 melt
#' @importFrom utils write.table
#' @export
#' @examples
#' suppressMessages(require(GenomicRanges))
#' suppressMessages(require(BiocFileCache))
#' 
#' # read CNV data
#' phenoFile <- system.file("extdata","AGP1_CNV.txt",package="netDx")
#' pheno   <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)
#' colnames(pheno)[1] <- "ID"
#' pheno <- pheno[!duplicated(pheno$ID),]
#' 
#' # create GRanges with patient CNVs
#' cnv_GR    <- GRanges(pheno$seqnames,IRanges(pheno$start,pheno$end),
#'                         ID=pheno$ID,LOCUS_NAMES=pheno$Gene_symbols)
#' 
#' # get gene coordinates
#' geneURL <- paste("http://download.baderlab.org/netDx/",
#' 	"supporting_data/refGene.hg18.bed",sep="")
#' cache <- rappdirs::user_cache_dir(appname = "netDx")
#' bfc <- BiocFileCache::BiocFileCache(cache,ask=FALSE)
#' rid_rec <- bfcquery(bfc, "hg18_genes", "rname")
#' rid <- rid_rec$rid
#' if (!length(rid)) {
#' 	rid <- names(bfcadd(bfc, "hg18_genes", geneURL))
#' }
#' if (!isFALSE(bfcneedsupdate(bfc, rid))){
#' 	bfcdownload(bfc, rid,ask=FALSE)
#' }
#' geneFile <- bfcrpath(bfc,rids=rid)
#' genes <- read.delim(geneFile,sep="\t",header=FALSE,as.is=TRUE)
#' genes <- genes[which(genes[,4]!=""),]
#' gene_GR     <- GRanges(genes[,1],IRanges(genes[,2],genes[,3]),
#'    name=genes[,4])
#' 
#' # create GRangesList of pathways
#' pathFile <- fetchPathwayDefinitions("February",2018,verbose=TRUE)
#' pathwayList <- readPathways(pathFile)
#' path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
#' 
#' #### uncomment to run - takes 5 min
#' #out <- buildPredictor_sparseGenetic(pheno, cnv_GR, "case",
#' #                             path_GRList,outDir,
#' #                             numSplits=3L, featScoreMax=3L,
#' #                             enrichLabels=TRUE,numPermsEnrich=20L,
#' #                             numCores=1L)
#' #summary(out)
#' #head(out$cumulativeFeatScores)
#' 
buildPredictor_sparseGenetic <- function(phenoDF,cnv_GR,predClass,
	group_GRList,outDir=tempdir(),
	numSplits=3L, featScoreMax=10L,
	filter_WtSum=100L,
	enrichLabels=TRUE,enrichPthresh=0.07,numPermsEnrich=2500L,minEnr=-1,
	numCores=1L,FS_numCores=NULL,...) {

	netDir <- sprintf("%s/networks_orig",outDir)

#message("making rangesets")
	netList <- makePSN_RangeSets(cnv_GR, group_GRList,netDir,
		verbose=FALSE)

#message("counting patients in net")
	p 	<- countPatientsInNet(netDir,netList, phenoDF$ID)
#message("updating nets")
	tmp	<- updateNets(p,phenoDF,writeNewNets=FALSE,verbose=FALSE)

	netmat	<- tmp[[1]]
	phenoDF	<- tmp[[2]] 

	if (is.null(FS_numCores)) FS_numCores <- max(1,numCores-1)
	# split into testing and training - resampling mode
	message("* Resampling train/test samples")
	TT_STATUS 	<- splitTestTrain_resampling(phenoDF, nFold=numSplits,
		predClass=predClass, verbose=TRUE)
	p_full 		<- netmat
	pheno_full	<- phenoDF
	if (any(colnames(pheno_full)%in% "TT_STATUS")) {
		message(paste("** Warning, found TT_STATUS column. ",
			"netDx adds its own column so this one will be removed **",
			sep=""))
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
							oldNetDir=netDir, newNetDir=trainNetDir,
							verbose=FALSE)
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
			sep=.Platform$file.sep)
		pathwayRank	<- compileFeatureScores(nrankFiles,
			filter_WtSum=filter_WtSum,verbose=FALSE)
		write.table(pathwayRank,file=sprintf("%s/pathwayScore.txt",resDir),
			col.names=TRUE,row.names=FALSE,quote=FALSE)
		
		pScore[[k]]	<- pathwayRank
	}
	phenoDF$TT_STATUS <- NULL

	# Measure performance after adding pathway tally across pplits
	out2 <- RR_featureTally(netmat, phenoDF, TT_STATUS, predClass,pScore,
			outDir,enrichLabels, enrichedNets,verbose=FALSE,
			maxScore=numSplits*featScoreMax)

	colnames(netmat) <- sub("_cont.txt","",colnames(netmat))
	
	out <- list(
		netmat=netmat,
		pheno=phenoDF,
		TT_STATUS=TT_STATUS,
		pathwayScores=pScore,
		enrichedNets=enrichedNets)	
	for (k in names(out2)) {
		out[[k]] <- out2[[k]]
	}

return(out)
}
