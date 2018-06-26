#' write pathways for selected nets to view as an EnrichmentMap in
#' Cytoscape
#'
require(netDx)
require(netDx.examples)
cutoff <- 9

datDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_170208"

# limit to genes in dataset
mutFile <- "/mnt/data2/BaderLab/PanCancer_GBM/input/from_firehose/GBM_core_somatic_mutations.txt"
xprFile <- "/mnt/data2/BaderLab/PanCancer_GBM/input/GBM_mRNA_core.txt"

xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]

for (gp in c("SURVIVENO","SURVIVEYES")) {
	cat(sprintf("%s\n",gp))
	scoreFile <- sprintf("%s/rng1/%s/GM_results/%s_pathway_CV_score.txt",
				 datDir,gp,gp)
	netScores <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
	           path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)
	
	netScores <- subset(netScores, netScores[,2]>=cutoff)

	outFile <- sprintf("GBM_somMut_%s_cutoff%i.gmt",basename(scoreFile),
					   cutoff)
	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))
	
	# RNA data
	idx <- grep(".profile",netScores[,1])
	netScores[,1] <- sub(".profile","",netScores[,1])	
	pList <- pathwayList[which(names(pathwayList)%in% netScores[,1])]	
	for (k in names(pList))  {
		cur <- pList[[k]]
		cur <- intersect(cur, xpr_genes) # limit to interrogated genes
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(cur,collapse="\t")),
			file=outFile,append=TRUE)
	}

	# somatic mutation pathways
	netScores <- netScores[-idx,]  # remove RNA variables.
	tmp <- sub("MUT_", "", netScores[,1])	
	tmp <- sub("_cont","",tmp)
	pList <- pathwayList[which(names(pathwayList) %in% tmp)]
	names(pList) <- paste("MUT",names(pList),sep="_")
	for (k in names(pList)) {
		cur <- pList[[k]]
		cur <- intersect(cur, mut_genes) # limit to interrogated genes
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(cur,collapse="\t")),
			file=outFile,append=TRUE)
	}
}
