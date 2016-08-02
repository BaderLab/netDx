#' write pathways for selected nets to view as an EnrichmentMap in
#' Cytoscape
#'
scoreFile <- "/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test/LumA_pathwayScore.txt"
netScores <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
outDir <- "/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test"
cutoff <- 13

require(netDx)
require(netDx.examples)
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

netScores <- subset(netScores, netScores[,2]>=cutoff)
netScores[,1] <- sub(".profile","",netScores[,1])

pList <- pathwayList[which(names(pathwayList)%in% netScores[,1])]

outFile <- sprintf("%s/%s_cutoff%i.gmt",outDir,basename(scoreFile),cutoff)
system(sprintf("touch %s",outFile))
for (k in names(pList)) 
	cat(sprintf("%s\t%s\t%s\n", k,k,paste(pList[[k]],collapse="\t")),
		file=outFile,append=TRUE)

