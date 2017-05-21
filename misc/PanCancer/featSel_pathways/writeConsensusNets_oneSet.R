# write names of nets that pass consensus across multiple runs
require(netDx)

#' @param datDir (char) path to output results directory. Directory structure
#' should be:  <datDir>/rng<runNumber>/<className>/GM_results/*_pathway_CV_score.txt
#' @param predClasses (char) classes that are being predicted
#' @param consCutoff (integer) nets must pass this cutoff in at least 
#' pctPass of the iterations to be written to file
#' @param pctPass (numeric between 0 and 1) fraction of iterations that
#' a net's score must pass consCutoff for, to be included in the consensus
#' map
#' @param netInfo (list)
#' keys are user-defined network types (e.g. clinical, mutations, expression),
#' value is a list with the entry:
#' 1) pattern: (char) a unique pattern to grep networks of this type in 
#' pathway_CV_score.txt. See note with ***.
#' 2) netSet: (list or char) key-value pair of network names and included units. 
#' e.g. pathway names and member genes.
#' If set to "base", the network name is written as is to the output file
#' *** Order in this list matters. If a net category does not have a unique
#' pattern, list it after others that do. Then those nets will be first
#' accounted for, and the residuals will be grep-ed for with this pattern.
#' For example, if mutation nets are MUT_*_cont and clinical nets have *_cont
#' then list the mutation ones first. After these are removed, any with *_cont
#' remaining should be clinical.
#' @param outPfx (char) prefix to output filename (absolute path)
#' @param trimFromName (char) strings to trim from net names before writing
#' to file
#' @param gmtFile (char) path to gmt file mapping set names to units
#' @param unitList (char) units measured in this predictor
#' @return no values. 
#' Writes several files with consensus/all net info
#' Outdir is <outPfx>_<gp>_thresh<i>_pctPass<1.2f>
#' 1) outDir/AllNets.txt  : All input nets
#' 2) outDir/consNets.txt : Consensus nets for cutoff
#' 3) outDir/netScores.txt : Table with net-level score
#' for all iterations provided
#' 4) outDir/consensus.gmt
#' 5) outDir/pheno.txt
#' 6) outDir/PSN.txt
#' 7) outDir/netInfo.txt
writeConsensusNets <- function(datDir,
	predClasses=c("SURVIVEYES","SURVIVENO"),netInfo,
	consCutoff=7L,pctPass=0.5,outPfx="./pred_",
	gmtFile,unitList) {

	outPfx <- sprintf("%s_thresh%i_pctPass%1.2f", outPfx,consCutoff, pctPass)
	for (gp in predClasses) {
		rngDirs <- dir(path=datDir, pattern="^rng")
		rngDirs <- setdiff(rngDirs, rngDirs[grep("tar.gz",rngDirs)])

	cat(sprintf("Got %i iterations\n", length(rngDirs)))
	netColl <- list()
	for (curDir in rngDirs) {
		scoreFile <- sprintf("%s/%s/%s/GM_results/%s_pathway_CV_score.txt",
				 datDir,curDir,gp,gp)
		tmp	 <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
		colnames(tmp)[1] <- "PATHWAY_NAME"
		netColl[[curDir]] <- tmp
	}
	# filter for nets meeting cutoff criteria
	cat("* Computing consensus\n")
	cons <- getNetConsensus(netColl); x1 <- nrow(cons)

	outFile <- sprintf("%s_%s_AllNets.txt",outPfx,gp)
	write.table(cons[,1],file=outFile,sep="\t",col=T,row=F,quote=F)

	outFile <- sprintf("%s_%s_netScores.txt",outPfx,gp)
	write.table(cons,file=outFile,sep="\t",col=T,row=F,quote=F)
	
	na_sum <- rowSums(is.na(cons))
	full_cons <- cons
	cons <- cons[which(na_sum < 1),]

	### uncomment this section to plot na_sum vs num_passes_cutoff
	 blah <- rowSums(full_cons[,-1]>=consCutoff,na.rm=T)
###	 plot(na_sum,blah, xlab="# splits for which net wasn't chosen", 
###			ylab="num splits for which net passes cutoff")
###	idx <- which(na_sum<1)
###	points(na_sum[idx],blah[idx],col='red')
###	abline(h=40,lty=3,col='red')
###	abline(h=50,lty=3,col='red')
###	title(basename(outPfx))
###	cat(sprintf("\t%i of %i scored in all rounds\n",nrow(cons),x1))

	avg_score <- rowMeans(cons[,-1])
	passes_cutoff <- cons[,-1]>= consCutoff;
	idx <- which(rowSums(passes_cutoff)>=round(pctPass*length(rngDirs)))
	cat(sprintf("\t** Consensus: %i nets scored >= %i **\n",
		length(idx),consCutoff))
	
	outFile <- sprintf("%s_%s_consNets.txt", outPfx, gp)
	cons <- cons[idx,]
	write.table(cons[,1],file=outFile,sep="\t",col=T,row=F,quote=F)
	}

	return(outPfx)
}
