# --------------------------------------------------------------
# writes high-scoring nets for each cancer type
rm(list=ls())
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_featSel_pathway_170426",rootDir)
)

source("writeConsensusNets_oneSet.R")
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/2017_PanCancer_Survival/featSel_pathways_170426/featSelNets", 
		rootDir)
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

logFile <- sprintf("%s/getConsNets_pathways_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
	for (curSet in "KIRC") { # names(dirList)){
		cat("-----------------------------\n")
		cat(sprintf("%s\n",curSet))
		cat("-----------------------------\n")
		for (cutoff in 9:10) {
		cat("-----------------------------\n")
		cat(sprintf("\tcutoff = %i\n", cutoff))
		cat("-----------------------------\n")

		outPfx<-sprintf("%s/%s_thresh%i",outDir,curSet,cutoff)
		writeConsensusNets(datDir=dirList[[curSet]],
			outPfx=outPfx,consCutoff=cutoff,pctPass=.75)

		cat("YES\n")
		dat <- read.delim(sprintf("%s_SURVIVEYES_consNets.txt",outPfx),
			sep="\t",header=T,as.is=T)
		print(dat)
		cat("\n")
		cat("NO\n")
		dat <- read.delim(sprintf("%s_SURVIVENO_consNets.txt",outPfx),
			sep="\t",header=T,as.is=T)
		print(dat)

		cat("-----------------------------\n")
		}
	}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
