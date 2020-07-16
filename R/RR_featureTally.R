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
#' @param enrichLabels (logical) was network label enrichment used?
#' @param enrichedNets (list of chars) networks passing network label 
#' enrichment
#' @param maxScore (integer) max achievable score for pathways
#' corresponding to N-way resampling
#' @param verbose (logical) print messages
#' @return (list) 
#' 1) cumulativeFeatScores: pathway name, cumulative score over 
#' N-way data resampling.
#' 2) performance_denAllNets: positive,negative calls at each cutoff: 
#' network score cutoff (score); num networks at cutoff (numPathways) ; 
#' total +, ground truth (pred_tot); + calls (pred_ol); 
#' + calls as pct of total (pred_pct); total -, ground truth (other_tot) ; 
#' - calls  (other_ol) ; - calls as pct of total (other_pct) ; ratio of 
#' pred_pct and other_pct (rr) ; min. pred_pct in all resamplings 
#' (pred_pct_min) ; max pred_pct in all resamplings (pred_pct_max) ; 
#' min other_pct in all resamplings (other_pct_min); max other_pct in all
#' resamplings (other_pct_max)
#' 3) performance_denEnrichedNets: positive, negative calls at each cutoff 
#' label enrichment option: format same as performance_denAllNets. 
#' However, the denominator here is limited to
#' patients present in networks that pass label enrichment
#' 4) resamplingPerformance: breakdown of performance for each of the 
#' resamplings, at each of the cutoffs. 
#' This is a list of length 2, one for allNets and one for enrichedNets. 
#' The value is a matrix with (resamp * 7) columns and S rows,
#' one row per score. The columns contain the following information
#' per resampling:
#' 1) pred_total: total num patients of predClass
#' 2) pred_OL: num of pred_total with a CNV in the selected net
#' 3) pred_OL_pct: 2) divided by 1) (percent)
#' 4) other_total: total num patients of other class(non-predClass)
#' 5) other_OL: num of other_total with CNV in selected net
#' 6) other_OL_pct: 5) divided by 4) (percent)
#' 7) relEnr: 6) divided by 3).
#' @importFrom utils write.table
#' @export 
#' @examples
#' data(cnv_patientNetCount) # patient presence/absence in nets
#' data(cnv_pheno)		# patient ID, label
#' data(cnv_netScores)	# network scores for resampling
#' data(cnv_TTstatus)	# train/test status
#' data(cnv_netPass) 	# nets passing label enrichment
#'
#' d <- tempdir()
#' out <- RR_featureTally(cnv_patientNetCount,
#' 		cnv_pheno,cnv_TTstatus,"case",cnv_netScores,
#' 		outDir=d,enrichLabels=TRUE,enrichedNets=cnv_netPass,
#' 		maxScore=30L)
#' print(summary(out))
RR_featureTally <- function(netmat,phenoDF,TT_STATUS,predClass,
	pScore,outDir=tempdir(),enrichLabels=TRUE,
	enrichedNets,maxScore=30L,
	verbose=FALSE) {

# tally pathway score across resamplings
pTally 			<- list()
for (k in seq_len(length(TT_STATUS))) {
	dat <- pScore[[k]]
	dat[,1] <- as.character(dat[,1])
	for (m in seq_len(nrow(dat))) {
		curp <- dat[m,1]
		if (!curp %in% names(pTally)) pTally[[curp]] <- 0
		pTally[[curp]] <- pTally[[curp]] + dat[m,2]
	}
}

# write pathway tally to file
pathDF <- data.frame(PATHWAY_NAME=names(pTally),SCORE=unlist(pTally))
pathDF[,2] <- as.integer(as.character(pathDF[,2]))
pathDF <- pathDF[order(pathDF[,2],decreasing=TRUE),]

tmpOut <- paste(outDir,"pathway_cumTally.txt",sep=getFileSep())
write.table(pathDF, file=tmpOut,sep="\t",col.names=TRUE, 
	row.names=FALSE,quote=FALSE)

out <- list()
tmp <- pathDF
tmp[,1] <- sub("_cont.txt","",tmp[,1])
out[["cumulativeFeatScores"]] <- tmp

# now run test
scoreColl <- seq_len(maxScore)
# two sets of results: one for den=allnets and one for den=enrichedNets
outdf <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
colnames(outdf) <- c("score", "numPathways",
	"pred_tot","pred_ol","pred_pct","other_tot","other_ol","other_pct", "rr",
	"pred_pct_min","pred_pct_max","other_pct_min","other_pct_max")

# and two more matrices for training samples
outdf_train <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
colnames(outdf_train) <- colnames(outdf)

if (enrichLabels) {
	outdf_enriched <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
	colnames(outdf_enriched) <- colnames(outdf)
	outdf_enriched_tr <- matrix(NA, nrow=length(scoreColl),ncol=9+4)
	colnames(outdf_enriched_tr) <- colnames(outdf)
}

ctr <- 1
# vector of patients contributing to OR. Union over all test resamplings
# order corresponds to that of the scoreColl vector of scores
predContr		<- rep("",length(scoreColl))
otherContr		<- rep("",length(scoreColl))
predContr_cl	<- rep("",length(scoreColl))
otherContr_cl	<- rep("",length(scoreColl))
resampPerf <- list(allNets=list(),enrichedNets=list())

for (setScore in scoreColl){
	selPath <- pathDF[which(pathDF[,2]>=setScore),1]
	## uncomment to test with original pathways
	if (verbose) message(sprintf("Thresh = %i ; %i pathways",
		setScore,length(selPath)))

	currmat				<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	currmat_enriched 		<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	# train
	currmat_train 		<- matrix(NA, nrow=length(TT_STATUS),ncol=7)
	currmat_enriched_tr 	<- matrix(NA, nrow=length(TT_STATUS),ncol=7)

	# store highest 
	predCurr	<- ""
	otherCurr	<- ""
	predCurr_cl <- ""
	otherCurr_cl <- ""

	for (k in seq_len(length(TT_STATUS))) {
		if (verbose) message(sprintf("\t(k = %i)",k))
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

		if (enrichLabels){
		# --------------------------------------------------------------
		# then run for denominator = label enriched nets
		pheno_test	<- phenoDF[which(TT_STATUS[[k]]%in% "TEST"),]
		p_test      <- netmat[which(rownames(netmat)%in% pheno_test$ID),]
		p_test		<- p_test[,
					which(colnames(p_test) %in% enrichedNets[[k]])]
		tmp 	<- updateNets(p_test,pheno_test,writeNewNets=FALSE,
							  verbose=FALSE)
		p_test	<- tmp[[1]]; pheno_test <- tmp[[2]];
		tmp <- getOR(p_test, pheno_test,predClass,selPath,verbose=FALSE)
		x <- tmp$stats
		# case: total, OL, pctOL. control: total, OL, pctOL. RelEnr.
		currmat_enriched[k,] <- c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],
					tmp$relEnr)

		predCurr_cl		<- c(predCurr_cl, intersect(tmp$OLsamps, 
			pheno_test$ID[which(pheno_test$STATUS %in% predClass)]))
		otherCurr_cl	<- c(otherCurr_cl, intersect(tmp$OLsamps,
			pheno_test$ID[which(!pheno_test$STATUS %in% predClass)]))

		# now for training samples
		pheno_train <- phenoDF[which(TT_STATUS[[k]]%in% "TRAIN"),]
		p_train     <- netmat[which(rownames(netmat)%in% pheno_train$ID),]
		p_train		<- p_train[,
					which(colnames(p_train) %in% enrichedNets[[k]])]
		tmp 	<- updateNets(p_train,pheno_train,writeNewNets=FALSE,
							  verbose=FALSE)
		p_train	<- tmp[[1]]; pheno_train <- tmp[[2]];
		tmp <- getOR(p_train, pheno_train,predClass,selPath,
					 verbose=FALSE)
		
		x <- tmp$stats
		currmat_enriched_tr[k,] <- 
			c(x[1,1],x[1,2],x[1,3],x[2,1],x[2,2],x[2,3],tmp$relEnr)
		}

	} # end loop over k resamplings for a given score

	# store resampling-wise info before averaging
	resampPerf[["allNets"]][[setScore]] <- currmat
	resampPerf[["enrichedNets"]][[setScore]] <- currmat_enriched

	predCurr 	<- unique(predCurr)
	otherCurr	<- unique(otherCurr) 
	if (verbose) {message(sprintf("\t# contrib: %i pred ; %i other",
				length(predCurr),length(otherCurr)))
	}

	predContr[ctr]	<- paste(predCurr,collapse=",")
	otherContr[ctr] <- paste(otherCurr,collapse=",")

	if (enrichLabels) {
		predCurr_cl 	<- unique(predCurr_cl)
		otherCurr_cl	<- unique(otherCurr_cl) 
		if (verbose) {
			message(paste("\tLABEL ENRICHMENT: ",
				sprintf("# contributing: %i pred ; %i other",
			length(predCurr_cl),length(otherCurr_cl)),sep=""))
		}
	
		predContr_cl[ctr]	<- paste(predCurr_cl,collapse=",")
		otherContr_cl[ctr]	<- paste(otherCurr_cl,collapse=",")
	}
	
	outdf[ctr,] <- c(setScore,length(selPath),colMeans(currmat),
		min(currmat[,3]),max(currmat[,3]),min(currmat[,6]),max(currmat[,6]))

	outdf_train[ctr,] <- c(setScore,length(selPath),colMeans(currmat_train),
		min(currmat_train[,3]),max(currmat_train[,3]),
		min(currmat_train[,6]),max(currmat_train[,6]))

	if (enrichLabels){
	outdf_enriched[ctr,] <- c(setScore,length(selPath),
		colMeans(currmat_enriched),
		min(currmat_enriched[,3]),max(currmat_enriched[,3]),
		min(currmat_enriched[,6]),max(currmat_enriched[,6]))

	outdf_enriched_tr[ctr,] <- c(setScore,length(selPath),
		colMeans(currmat_enriched_tr),
		min(currmat_enriched_tr[,3]),max(currmat_enriched_tr[,3]),
		min(currmat_enriched_tr[,6]),max(currmat_enriched_tr[,6]))
	}

	ctr <- ctr+1
} # end loop over score cutoffs

numresamp <- nrow(resampPerf[[1]][[1]])
for (k in seq_len(length(resampPerf))) {
	tmp <- lapply(resampPerf[[k]], function(x) {as.numeric(t(x))})
	tmp <- do.call(rbind,tmp)
	x <- rep(c("pred_total","pred_OL","pred_OL_pct",
				"other_total","other_OL","other_OL_pct","relEnr"),
			numresamp)
	colnames(tmp) <- paste(x,rep(seq_len(numresamp),each=7),sep="_")
	rownames(tmp) <- scoreColl
	resampPerf[[k]] <- tmp
}

out[["resamplingPerformance"]] <- resampPerf

save(resampPerf,file=paste(outDir,"resamplingPerf.Rdata",sep=getFileSep()))

# add IDs of contributing samples
outdf <- data.frame(outdf)
outdf <- cbind(outdf, CONTRIBUT_PRED=predContr,
			   CONTRIBUT_OTHER=otherContr)

outFile <- paste(outDir,"RR_changeNetSum_stats_denAllNets.txt",
	sep=getFileSep())
write.table(outdf,file=outFile,sep="\t",
	col.names=TRUE,row.names=FALSE,quote=FALSE)
out[["performance_denAllNets"]] <- outdf

write.table(outdf_train,file=outFile,sep="\t",
		col.names=TRUE,row.names=FALSE,quote=FALSE)
out[["performance_denAllNets_TrainingSamples"]] <- outdf_train

if (enrichLabels) {
	# add IDs of contributing samples
	outdf_enriched <- data.frame(outdf_enriched)
	outdf_enriched <- cbind(outdf_enriched, CONTRIBUT_PRED=predContr_cl,
			   CONTRIBUT_OTHER=otherContr_cl)
	
	outFile <- paste(outDir,"RR_changeNetSum_stats_denEnrichedNets.txt",
		sep=getFileSep())
write.table(outdf_enriched,file=outFile,sep="\t",
	col.names=TRUE,row.names=FALSE,quote=FALSE)
out[["performance_denEnrichedNets"]] <- outdf_enriched
}

return(out)
}
