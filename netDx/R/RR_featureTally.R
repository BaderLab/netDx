#' Computes positive and negative calls upon changing stringency of 
#' feature selected networks (binary networks only)
#'
#' @details This function computes predictor performance in the context
#' of binary networks, where + and - calls are based on membership (or 
#' lack thereof) in feature selected networks. An example would be networks
#' based on CNV occurrence in cellular pathways; in this use case, a +
#' is based on patient membership in feature-selected networks. This
#' function takes the output data from a feature selection exercise
#' and computes the number and fraction of positive and negative calls at
#' each level of feature selection stringency. The output of this function
#' can then be used to compute performance measures such as the ROC or 
#' precision-recall curve.
#' @param netmat (matrix) output of countPatientsInNet. Should contain
#' all patients in dataset that overlap 1+ network
#' @param phenoDF (data.frame) patient ID and STATUS
#' @param TT_STATUS (list) output of splitTestTrain_partition; should be
#' same as used for cross validation
#' @param predClass (char) class to be predicted
#' @param pScore (list of data.frames) contains 10-fold CV score, one entry
#' for each resampling of the data. The data.frame has two columns: 
#' 1) pathway name, 2) pathway score
#' @param outDir (char) path to dir where results should be written
#' @param cliqueFilter (logical) was clique filtering used?
#' @param cliqueNets (list of chars) networks passing clique filtering
#' @param maxScore (integer) max achievable score for pathways
#' corresponding to N-way resampling
#' @param verbose (logical) print messages
#' @return No value. Side effect of writing the following files: 
#' 1) pathway tally file: \code{<outDir>/pathway_cumTally.txt} : 
#' pathway name, cumulative score over N-way data resampling.
#' 2) positive,negative calls at each cutoff: 
#' \code{<outDir>/RR_changeNetSum_stats_denAllNets.txt}: 
#' network score cutoff (score); num networks at cutoff (numPathways) ; 
#' total +, ground truth (pred_tot); + calls (pred_ol); 
#' + calls as pct of total (pred_pct); total -, ground truth (other_tot) ; 
#' - calls  (other_ol) ; - calls as pct of total (other_pct) ; ratio of 
#' pred_pct and other_pct (rr) ; min. pred_pct in all resamplings 
#' (pred_pct_min) ; max pred_pct in all resamplings (pred_pct_max) ; 
#' min other_pct in all resamplings (other_pct_min); max other_pct in all
#' resamplings (other_pct_max)
#' 3) positive, negative calls at each cutoff ; clique filtering option:  
#' \code{<outDir>/RR_changeNetSum_stats_denCliqueNets.txt}: format same
#' as the previous file. However, the denominator here is limited to
#' patients present in networks that pass clique filtering. 
#' 4) breakdown of performance for each of the resamplings, at each of the 
#' cutoffs: <outDir>/resamplingPerf.Rdata: list of length 2, one for allNets and
#' one for cliqueNets. The value is a matrix with (resamp * 7) columns and S rows,
#' one row per score. The columns contain the followin information per resampling:
#' 1) pred_total: total num patients of predClass
#' 2) pred_OL: num of pred_total with a CNV in the selected net
#' 3) pred_OL_pct: 2) divided by 1) (percent)
#' 4) other_total: total num patients of other class(non-predClass)
#' 5) other_OL: num of other_total with CNV in selected net
#' 6) other_OL_pct: 5) divided by 4) (percent)
#' 7) relEnr: 6) divided by 3).
#' @importFrom utils write.table
#' @export 
RR_featureTally <- function(netmat,phenoDF,TT_STATUS,predClass,
	pScore,outDir,cliqueFilter=TRUE, cliqueNets,maxScore=30L,
	verbose=FALSE) {

# tally pathway score across resamplings
testMode <- FALSE # when true doesn't write files.

pTally 			<- list()
for (k in 1:length(TT_STATUS)) {
	dat <- pScore[[k]]
	dat[,1] <- as.character(dat[,1])
	for (m in 1:nrow(dat)) {
		curp <- dat[m,1]
		if (!curp %in% names(pTally)) pTally[[curp]] <- 0
		pTally[[curp]] <- pTally[[curp]] + dat[m,2]
	}
}

# write pathway tally to file
pathDF <- data.frame(PATHWAY_NAME=names(pTally),SCORE=unlist(pTally))
pathDF[,2] <- as.integer(as.character(pathDF[,2]))
pathDF[,1] <- paste(as.character(pathDF[,1]),".txt",sep="")
pathDF <- pathDF[order(pathDF[,2],decreasing=TRUE),]

tmpOut <- sprintf("%s/pathway_cumTally.txt",outDir)
if (!testMode){
write.table(pathDF, file=tmpOut,sep="\t",col=TRUE, row=FALSE,quote=FALSE)
}

# now run test
scoreColl <- 1:maxScore
# two sets of results: one for den=allnets and one for den=cliquenets
outdf <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
colnames(outdf) <- c("score", "numPathways",
	"pred_tot","pred_ol","pred_pct","other_tot","other_ol","other_pct", "rr",
	"pred_pct_min","pred_pct_max","other_pct_min","other_pct_max")

# and two more matrices for training samples
outdf_train <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
colnames(outdf_train) <- colnames(outdf)

if (cliqueFilter) {
	outdf_clique <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
	colnames(outdf_clique) <- colnames(outdf)
	outdf_clique_tr <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
	colnames(outdf_clique_tr) <- colnames(outdf)
}

ctr <- 1
# vector of patients contributing to OR. Union over all test resamplings
# order corresponds to that of the scoreColl vector of scores
predContr		<- rep("",length(scoreColl))
otherContr		<- rep("",length(scoreColl))
predContr_cl	<- rep("",length(scoreColl))
otherContr_cl	<- rep("",length(scoreColl))
resampPerf <- list(allNets=list(),cliqueNets=list())

for (setScore in scoreColl){
	selPath <- pathDF[which(pathDF[,2]>=setScore),1]
	## uncomment to test with original pathways
	if (verbose) cat(sprintf("Thresh = %i ; %i pathways\n",setScore,length(selPath)))

	currmat				<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	currmat_clique 		<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	# train
	currmat_train 		<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	currmat_clique_tr 	<- matrix(NA, nrow=length(TT_STATUS),ncol=7)

	# store highest 
	predCurr	<- ""
	otherCurr	<- ""
	predCurr_cl <- ""
	otherCurr_cl <- ""

	for (k in 1:length(TT_STATUS)) {
		if (verbose) cat(sprintf("\t(k = %i)",k))
		# --------------------------------------------------------------
		# first run for denominator = all nets
		# set pheno and p to contain only test samples
		pheno_test	<- phenoDF[which(TT_STATUS[[k]]%in% "TEST"),]
		p_test      <- netmat[which(rownames(netmat)%in% pheno_test$ID),]
		tmp 		<- updateNets(p_test,pheno_test,writeNewNets=FALSE,
							  verbose=FALSE)
		p_test	<- tmp[[1]]; pheno_test <- tmp[[2]];

		tmp <- getOR(p_test, pheno_test,predClass,selPath,verbose=FALSE)

		predCurr	<- c(predCurr, intersect(tmp$OLsamps, 
			pheno_test$ID[which(pheno_test$STATUS %in% predClass)]))
		otherCurr	<- c(otherCurr, intersect(tmp$OLsamps,
			pheno_test$ID[which(!pheno_test$STATUS %in% predClass)]))
		
		x <- tmp$stats
		currmat[k,] <- c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],
						tmp$relEnr)

		rm(x,tmp,pheno_test,p_test)

		# now training samples
		pheno_train	<- phenoDF[which(TT_STATUS[[k]]%in% "TRAIN"),]
		p_train     <- netmat[which(rownames(netmat)%in% pheno_train$ID),]
		tmp 		<- updateNets(p_train,pheno_train,writeNewNets=FALSE,
							  verbose=FALSE)
		p_train	<- tmp[[1]]; pheno_train <- tmp[[2]];

		tmp <- getOR(p_train, pheno_train,predClass,selPath,
					 verbose=FALSE)
		x <- tmp$stats
		currmat_train[k,] <- c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],
						tmp$relEnr)
		rm(x,tmp,pheno_train,p_train)

		if (cliqueFilter){
		# --------------------------------------------------------------
		# then run for denominator = clique filtered nets
		pheno_test	<- phenoDF[which(TT_STATUS[[k]]%in% "TEST"),]
		p_test      <- netmat[which(rownames(netmat)%in% pheno_test$ID),]
		p_test		<- p_test[,
					which(colnames(p_test) %in% cliqueNets[[k]])]
		tmp 	<- updateNets(p_test,pheno_test,writeNewNets=FALSE,
							  verbose=FALSE)
		p_test	<- tmp[[1]]; pheno_test <- tmp[[2]];
		tmp <- getOR(p_test, pheno_test,predClass,selPath,verbose=FALSE)
		x <- tmp$stats
		# case: total, OL, pctOL. control: total, OL, pctOL. RelEnr.
		currmat_clique[k,] <- c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],
					tmp$relEnr)

		predCurr_cl		<- c(predCurr_cl, intersect(tmp$OLsamps, 
			pheno_test$ID[which(pheno_test$STATUS %in% predClass)]))
		otherCurr_cl	<- c(otherCurr_cl, intersect(tmp$OLsamps,
			pheno_test$ID[which(!pheno_test$STATUS %in% predClass)]))

		# now for training samples
		pheno_train <- phenoDF[which(TT_STATUS[[k]]%in% "TRAIN"),]
		p_train     <- netmat[which(rownames(netmat)%in% pheno_train$ID),]
		p_train		<- p_train[,
					which(colnames(p_train) %in% cliqueNets[[k]])]
		tmp 	<- updateNets(p_train,pheno_train,writeNewNets=FALSE,
							  verbose=FALSE)
		p_train	<- tmp[[1]]; pheno_train <- tmp[[2]];
		tmp <- getOR(p_train, pheno_train,predClass,selPath,
					 verbose=FALSE)
		
		x <- tmp$stats
		currmat_clique_tr[k,] <- 
			c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],tmp$relEnr)
		}

	} # end loop over k resamplings for a given score

	# store resampling-wise info before averaging
	resampPerf[["allNets"]][[setScore]] <- currmat
	resampPerf[["cliqueNets"]][[setScore]] <- currmat_clique

	predCurr 	<- unique(predCurr)
	otherCurr	<- unique(otherCurr) 
	if (verbose) {cat(sprintf("\t# contrib: %i pred ; %i other\n",
				length(predCurr),length(otherCurr)))
	}

	predContr[ctr]	<- paste(predCurr,collapse=",")
	otherContr[ctr] <- paste(otherCurr,collapse=",")

	if (cliqueFilter) {
		predCurr_cl 	<- unique(predCurr_cl)
		otherCurr_cl	<- unique(otherCurr_cl) 
		if (verbose) {
			cat(sprintf("\tCLIQUE: # contributing: %i pred ; %i other\n",
			length(predCurr_cl),length(otherCurr_cl)))
		}
	
		predContr_cl[ctr]	<- paste(predCurr_cl,collapse=",")
		otherContr_cl[ctr]	<- paste(otherCurr_cl,collapse=",")
	}
	
	if (verbose) cat("\n")
	outdf[ctr,] <- c(setScore,length(selPath),colMeans(currmat),
		min(currmat[,3]),max(currmat[,3]),min(currmat[,6]),max(currmat[,6]))

	outdf_train[ctr,] <- c(setScore,length(selPath),colMeans(currmat_train),
		min(currmat_train[,3]),max(currmat_train[,3]),
		min(currmat_train[,6]),max(currmat_train[,6]))

	if (cliqueFilter){
	outdf_clique[ctr,] <- c(setScore,length(selPath),
		colMeans(currmat_clique),
		min(currmat_clique[,3]),max(currmat_clique[,3]),
		min(currmat_clique[,6]),max(currmat_clique[,6]))

	outdf_clique_tr[ctr,] <- c(setScore,length(selPath),
		colMeans(currmat_clique_tr),
		min(currmat_clique_tr[,3]),max(currmat_clique_tr[,3]),
		min(currmat_clique_tr[,6]),max(currmat_clique_tr[,6]))
	}


	ctr <- ctr+1
	if (verbose) cat("\n")
} # end loop over score cutoffs

numresamp <- nrow(resampPerf[[1]][[1]])
for (k in 1:length(resampPerf)) {
	tmp <- lapply(resampPerf[[k]], function(x) {as.numeric(t(x))})
	tmp <- do.call(rbind,tmp)
	x <- rep(c("pred_total","pred_OL","pred_OL_pct",
				"other_total","other_OL","other_OL_pct","relEnr"),
			numresamp)
	colnames(tmp) <- paste(x,rep(1:numresamp,each=7),sep="_")
	rownames(tmp) <- scoreColl
	resampPerf[[k]] <- tmp
}

save(resampPerf,file=sprintf("%s/resamplingPerf.Rdata",outDir))

# add IDs of contributing samples
outdf <- data.frame(outdf)
outdf <- cbind(outdf, CONTRIBUT_PRED=predContr,
			   CONTRIBUT_OTHER=otherContr)

outFile <- sprintf("%s/RR_changeNetSum_stats_denAllNets.txt",
				   outDir)
if (!testMode){
write.table(outdf,file=outFile,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}

outFile <- sprintf("%s/RR_changeNetSum_stats_denAllNets_train.txt",
				   outDir)
if (!testMode){
write.table(outdf_train,file=outFile,sep="\t",
			col=TRUE,row=FALSE,quote=FALSE)
}

if (cliqueFilter) {
	# add IDs of contributing samples
	outdf_clique <- data.frame(outdf_clique)
	outdf_clique <- cbind(outdf_clique, CONTRIBUT_PRED=predContr_cl,
			   CONTRIBUT_OTHER=otherContr_cl)
	
	outFile <- sprintf("%s/RR_changeNetSum_stats_denCliqueNets.txt",
				   outDir)
	if (!testMode){
	write.table(outdf_clique,file=outFile,sep="\t",
		col=TRUE,row=FALSE,quote=FALSE)
	}

	outFile <- sprintf("%s/RR_changeNetSum_stats_denCliqueNets_train.txt",
				   outDir)
	if (!testMode){
	write.table(outdf_clique_tr,file=outFile,sep="\t",
		col=TRUE,row=FALSE,quote=FALSE)
	}
}
}
