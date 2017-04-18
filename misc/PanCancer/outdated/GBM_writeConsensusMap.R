#' write an enrichment map of the nets that on average score above
#' the cutoff, across a set of iterations.
#'
require(netDx)
require(netDx.examples)
cutoff <- 7
numRuns <- 35

datDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_noMut_170208"

# setup pathway info
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# collect nets for each group
for (gp in c("SURVIVENO","SURVIVEYES")) {
	cur <- list()
	for (k in 1:numRuns) {
		scoreFile <- sprintf("%s/rng%i/%s/GM_results/%s_pathway_CV_score.txt",
				 datDir,k,gp,gp)
		tmp	 <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
		colnames(tmp)[1] <- "PATHWAY_NAME"
		cur[[sprintf("rng%i",k)]] <- tmp
	}
	
	# keep only nets that were scored all the time and that
	# had a median score of 8.
	cons <- getNetConsensus(cur); x1 <- nrow(cons)
	na_sum <- rowSums(is.na(cons))
	cons <- cons[which(na_sum < 1),]
	cat(sprintf("\t%i of %i scored in all rounds\n",length(cons),x1))
	
	avg_score <- rowMeans(cons[,-1])
	passes_cutoff <- cons[,-1]>= cutoff;
	idx <- which(rowSums(passes_cutoff)>=round(numRuns/2))
	cat(sprintf("%i nets scored >= %i\n",length(idx),cutoff))
	print(cons[idx,1])
	cons <- cons[idx,]
	cons[,1] <- sub(".profile","",cons[,1])

	outFile <- sprintf("%s/%s_cutoff%i.txt",datDir, basename(scoreFile),
		cutoff)
	write.table(cons[,1],file=outFile,sep="\t",col=F,row=F,quote=F)


	pList <- pathwayList[which(names(pathwayList)%in% cons[,1])]
	outFile <- sprintf("%s/%s_cutoff%i.gmt",datDir,basename(scoreFile),
					   cutoff)
	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))
	for (k in names(pList)) 
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(pList[[k]],collapse="\t")),
		file=outFile,append=TRUE)
}

###for (gp in c("SURVIVENO","SURVIVEYES")) {
###	cat(sprintf("%s\n",gp))
###	
###	netScores <- subset(netScores, netScores[,2]>=cutoff)
###	netScores[,1] <- sub(".profile","",netScores[,1])
###	
###	pList <- pathwayList[which(names(pathwayList)%in% netScores[,1])]
###	
###}
