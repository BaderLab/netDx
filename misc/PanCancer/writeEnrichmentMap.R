#' write pathways for selected nets to view as an EnrichmentMap in
#' Cytoscape
#'
require(netDx)
require(netDx.examples)
cutoff <- 9

datDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_170208"

for (gp in c("SURVIVENO","SURVIVEYES")) {
	cat(sprintf("%s\n",gp))
	scoreFile <- sprintf("%s/%s/GM_results/%s_pathway_CV_score.txt",
				 datDir,gp,gp)
	netScores <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
	           path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)
	
	netScores <- subset(netScores, netScores[,2]>=cutoff)
	netScores[,1] <- sub(".profile","",netScores[,1])
	
	pList <- pathwayList[which(names(pathwayList)%in% netScores[,1])]
	
	outFile <- sprintf("%s/%s_cutoff%i.gmt",datDir,basename(scoreFile),
					   cutoff)
	if (file.exists(outFile)) unlink(outFile)
	system(sprintf("touch %s",outFile))
	for (k in names(pList)) 
		cat(sprintf("%s\t%s\t%s\n", k,k,paste(pList[[k]],collapse="\t")),
		file=outFile,append=TRUE)
}
