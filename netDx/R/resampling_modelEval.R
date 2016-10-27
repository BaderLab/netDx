#' Resampling-based predictor, evaluate test samples
#'
#' @param pheno (data.frame) patient metadata. ID, class labels (STATUS),
#' and whether train/test status (TT_STATUS): a patient is query if "TRAIN",
#' else is "TEST"
#' @param pdat (data.frame) patient data. Rows are unit variables, columns 
#' are patients. Rows and columns must be named.
#' @param predClass (char) class of interest to predict.
#' @param unitSets (list) groupings of unit variables; each grouping will
#' be converted into its own PSN. 
#' @param p_GR (GRanges) GRanges of patient CNVs. Has ID column in
#' metadata, containing patient IDs
#' @param unitSet_GR (list) sets of GRanges to group CNVs (e.g.
#' could have one GRanges per pathway, corresponding to regions in that 
#' pathway
#' @param bestCutoff (integer) nets with score >= bestCutoff are used in 
#' the final model
#' @param patNets (list of chars) for each class (key), subset of "unitSets"
#' for which patient nets must be created. e.g. feature-selected networks
#' for the corresponding patient subtype.
#' @param dataDir (char) dir containing subdirs corresponding to each 
#' subtype
#' @param outDir (char) directory in which results must be written
#' @param netScores (list) one key per subtype, and value is a data.frame
#' with net name (column 1) and score (column 2). If NULL, is obtained
#' using dataDir
#' @param numCores (integer) number of cores for parallel processing
#' @return (list)
#' 1. test_res: output of GM_predClass_once
#' 2. netScores (list) keys are subtypes, values are feature selected
#' nets
#' @export
resampling_modelEval <- function(pheno,pdat,predClass,unitSets,
	p_GR=NULL, unitSet_GR=NULL,bestCutoff,
	dataDir,outDir,netScores=NULL,numCores) {

	cat("* Evaluating model on test cases\n")
	subtypes <- unique(pheno$STATUS)
	
	if (is.null(netScores)) {
		cat("netScores not provided; extracting\n")
		##  pathways by subtype
		netScores <- list()
		for (g in subtypes) {
			pTally <- read.delim(
				sprintf("%s/%s_pathwayScore.txt",outDir,g),
				sep="\t",h=T,as.is=T)
			pTally <- pTally[which(pTally[,2]>=1),]
			pTally[,1] <- sub(".profile","",pTally[,1])
			pTally[,1] <- sub("_cont","",pTally[,1])
			netScores[[g]] <- pTally
		}
	}

	finalNets <- lapply(netScores, function(x) 
		x[which(x[,2]>=bestCutoff),1])

	testRes <- GM_predClass_once(pheno,pdat,predClass,unitSets,
	p_GR=p_GR, unitSet_GR=unitSet_GR,
	patNets=finalNets,outDir=outDir,numCores=numCores)

	return(list(testRes=testRes,netScores=netScores))	
}
