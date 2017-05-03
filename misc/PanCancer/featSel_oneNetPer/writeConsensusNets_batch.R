# --------------------------------------------------------------
# writes high-scoring nets for each cancer type
rm(list=ls())
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_oneNetPer_170426",rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/output/GBM_oneNetPer_170425",rootDir),
	OV=sprintf("%s/2017_TCGA_OV/output/OV_oneNetPer_170425",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/LUSC_oneNetPer_170425",rootDir)
)

source("writeConsensusNets_oneSet.R")
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/2017_PanCancer_Survival/oneNetPer_FeatSel/featSelNets", 
		rootDir)
if (!file.exists(outDir)) dir.create(outDir)

logFile <- sprintf("%s/getConsNets_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
	for (curSet in "KIRC") { # names(dirList)){
		cat("-----------------------------\n")
		cat(sprintf("%s\n",curSet))
		cat("-----------------------------\n")
		for (cutoff in 9:10) {
		cat(sprintf("\tcutoff = %i\n", cutoff))
		writeConsensusNets(datDir=dirList[[curSet]],
			outPfx=sprintf("%s/%s_thresh%i",outDir,curSet,cutoff),
			consCutoff=cutoff,pctPass=1)
		}
	}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
