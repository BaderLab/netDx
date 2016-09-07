#' write pathways for selected nets to view as an EnrichmentMap in
#' Cytoscape
#'
require(netDx)
require(netDx.examples)

## ### resampling result
scoreFile <- "/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test/LumA_pathwayScore.txt"
outDir <- "/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test"
cutoff <- 13

## ### one-time 10-fold CV
#scoreFile <- "/home/spai/tmp/TCGA_BRCA/LumA/GM_results/LumA_pathway_CV_score.txt"
#outDir <- "/home/spai/tmp/TCGA_BRCA"
#cutoff <- 9

netScores <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

netScores <- subset(netScores, netScores[,2]>=cutoff)
netScores[,1] <- sub(".profile","",netScores[,1])

pList <- pathwayList[which(names(pathwayList)%in% netScores[,1])]

outFile <- sprintf("%s/%s_cutoff%i.gmt",outDir,basename(scoreFile),cutoff)
if (file.exists(outFile)) unlink(outFile)
system(sprintf("touch %s",outFile))
for (k in names(pList)) 
	cat(sprintf("%s\t%s\t%s\n", k,k,paste(pList[[k]],collapse="\t")),
		file=outFile,append=TRUE)

idx <- grep("CNV_",netScores[,1])
if (any(idx)) {
	net <- netScores[idx,]
	net[,1] <- sub("_cont","",net[,1])
	net[,1] <- sub("CNV_","",net[,1])
	pList <- pathwayList[which(names(pathwayList)%in% net[,1])]
	for (k in names(pList)) 
		cat(sprintf("CNV_%s\tCNV_%s\t%s\n", k,k,
			paste(pList[[k]],collapse="\t")),
			file=outFile,append=TRUE)
}

