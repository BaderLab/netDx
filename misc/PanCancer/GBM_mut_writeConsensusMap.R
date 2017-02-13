#' write an enrichment map of the nets that on average score above
#' the cutoff, across a set of iterations.
#'
require(netDx)
require(netDx.examples)
cutoff <- 7
numRuns <- 30

datDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_170209"

# limit to genes in dataset
mutFile <- "/mnt/data2/BaderLab/PanCancer_GBM/input/from_firehose/GBM_core_somatic_mutations.txt"
xprFile <- "/mnt/data2/BaderLab/PanCancer_GBM/input/GBM_mRNA_core.txt"

xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]

# setup pathway info
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# collect nets for each group
for (gp in c("SURVIVENO","SURVIVEYES")) {
cat("----------------------\n")
cat(sprintf("%s\n",gp))
cat("----------------------\n")
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
	cat(sprintf("\t%i nets scored >= %i\n",length(idx),cutoff))
	#print(cons[idx,1])
	cons <- cons[idx,1]

	outFile <- sprintf("%s/%s_consensusNets_cutoff%i.txt",datDir,
		basename(scoreFile),cutoff)
	write.table(cons,file=outFile,sep="\t",col=F,row=F,quote=F)

	outFile <- sprintf("%s/%s_somMut_cutoff%i.gmt",datDir,
		basename(scoreFile),cutoff)

	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))

	# first write RNA-based pathways
	idx <- grep(".profile",cons)
	cons<- sub(".profile","",cons)
	pList <- pathwayList[which(names(pathwayList)%in% cons[idx])]
	cat(sprintf("\t > %i gene-expression based pathways\n",length(pList)))
	for (k in names(pList))  {
		cur <- pList[[k]]
		cur <- intersect(cur, xpr_genes) # limit to interrogated genes
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(cur,collapse="\t")),
			file=outFile,append=TRUE)
	}

	# somatic mutation pathways
	cons <- cons[-idx]  # remove RNA variables.
	tmp <- sub("MUT_", "", cons)
	tmp <- sub("_cont","",tmp)
	pList <- pathwayList[which(names(pathwayList) %in% tmp)]
	names(pList) <- paste("MUT",names(pList),sep="_")
	cat(sprintf("\t > %i mutation based pathways\n",length(pList)))
	for (k in names(pList)) {
		cur <- pList[[k]]
		cur <- intersect(cur, mut_genes) # limit to interrogated genes
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(cur,collapse="\t")),
			file=outFile,append=TRUE)
	}
	
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
