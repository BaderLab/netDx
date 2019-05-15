#' Evaluate model performance for various network cutoffs
#'
#' @details Currently only works when all data can be provided as a 
#' single matrix. Also currently works only for binary classification 
#' although extension to 3+-way classificaiton is likely straightforward.
#' @param pheno (data.frame) patient metadata. Patient ID (ID), subtypes
#' (STATUS), and whether they should be part of query or not (TT_STATUS; 
#' patients with TT_STATUS=train will be part of GM query)
#' @param pdat (matrix) patient data to build networks from. Rows are 
#' patients, columns are unit measures
#' @param p_GR (GRanges) GRanges of patient CNVs. Has ID column in
#' metadata, containing patient IDs. If NULL, assumes there is no
#' patient-range type data
#' @param unitSet_GR (list) sets of GRanges to group CNVs (e.g.
#' could have one GRanges per pathway, corresponding to regions in that 
#' pathway
#' @param netScores (list) scores of individual networks for each patient
#' label. Key is patient label; value is data.frame with two columns, 
#' PATHWAY_NAME and SCORE. PATHWAY_NAME should match names in unitSets
#' @param unitSets (list) unit groupings, each of which will be converted to
#' its own patient similarity network. 
#' @param maxScore (integer) maximum score achievable by a network
#' @param outDir (char) directory to store results in
#' @param numCores (integer) num cores for parallel processing
#' @param ... params for makePSN_NamedMatrix
#' @export
evalPerfOverCutoffs <- function(pheno,pdat,p_GR,unitSet_GR,
	predClass,netScores,unitSets,
	maxScore,outDir,numCores=1L,...) {

predRes <- list() 	## predRes[[cutoff]] contains predictions for 
					## score >= cutoff.
for (cutoff in 1:maxScore) predRes[[cutoff]] <- list()
subtypes <- names(netScores)
cat(sprintf("Got %i subtypes: { %s }\n", length(subtypes),
	paste(subtypes,collapse=",")))

for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g)
	dir.create(pDir)
	pTally <- netScores[[g]]
	pTally <- pTally[which(pTally[,2]>=1),]
	cat(sprintf("%s: %i pathways\n",g,nrow(pTally)))

	# create new db
	profDir <- sprintf("%s/profiles",pDir)
	tmp <- makePSN_NamedMatrix(pdat,rownames(pdat),
		unitSets[which(names(unitSets)%in% pTally[,1])],
		profDir,verbose=F,numCores=numCores,writeProfiles=TRUE ,...)

	if (!is.null(p_GR)) {
		cat("Got patient ranges - creating range-set nets\n")
		idx <- which(names(unitSet_GR) %in% sub("_cont","",pTally[,1]))
		if (length(idx)>0) {
			# note that range-sets with < 2 patients 
			# automatically don't get included
			netList2 <- makePSN_RangeSets(p_GR, 
				unitSet_GR[idx],profDir,verbose=FALSE)
		} else {
			cat("Not making any range-related nets\n")
		}
	}
	dbDir <- compileFeatures(profDir,pheno$ID,pDir,numCores=numCores)

	# query of all training samples for this class
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
			 pheno$TT_STATUS%in%"TRAIN")]

	# using new Genemania implementation
	for (cutoff in 1:maxScore){
	  cat(sprintf("\tCutoff = %i\n",cutoff))
	  qFile <- sprintf("%s/%s_cutoff%i",pDir,g,cutoff)
	  curr_p <- pTally[which(pTally[,2]>=cutoff),1]
	  idx <- grep("_cont$", curr_p,invert=TRUE)
	  if (length(idx)>0) {
	    curr_p[idx] <- paste(curr_p[idx],".profile",sep="")
	  }
	  if (length(curr_p)>0){
	    writeQueryFile(qSamps,curr_p,nrow(pheno),qFile)
	    resFile <- netDx::runQuery(dbDir$dbDir,qFile,resDir=pDir, numCores=numCores)
	    # system(sprintf("unlink %s", resFile)) # file does not exist
	  }
	}
	# collect rankings
	for (cutoff in 1:maxScore) {
		resFile <- sprintf("%s/%s_cutoff%i-results.report.txt.PRANK",
			pDir,g,cutoff)
		predRes[[cutoff]][[g]] <- getPatientRankings(resFile,pheno,g)
	}
}

# get rankings + compile performance at each cutoff
outClass <- list()
outmat <- matrix(NA, nrow=maxScore, ncol=9)
colnames(outmat) <- c("score","tp","tn","fp","fn","tpr","fpr",
					  "accuracy","ppv")

for (cutoff in 1:maxScore) {
	
	outClass[[cutoff]] <- predictPatientLabels(predRes[[cutoff]])
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

out <- list(predRes=predRes,predLabels=outClass,confmat=outmat)
save(out, file=sprintf("%s/predictionResults.Rdata", outDir))

out
}

