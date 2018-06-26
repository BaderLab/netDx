# write names of nets that pass consensus across multiple runs
require(netDx)

writeConsensusNets <- function(datDir,consCutoff=7L,
	predClasses=c("SURVIVEYES","SURVIVENO"),
	pctPass=0.5,outPfx="./pred_") {
	for (gp in predClasses) {
		rngDirs <- dir(path=datDir, pattern="^rng")
		rngDirs <- setdiff(rngDirs, rngDirs[grep("tar.gz",rngDirs)])

	cat(sprintf("Got %i iterations\n", length(rngDirs)))
	netColl <- list()
	for (curDir in rngDirs) {
		scoreFile <- sprintf("%s/%s/rna/%s/GM_results/%s_pathway_CV_score.txt",
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
}
