# --------------------------------------------------------------
# writes high-scoring nets for each cancer type
rm(list=ls())
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	#KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_featSel_pathway_170426",rootDir)
#	KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinNets_170430",rootDir)
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/pathway_170502",rootDir)
)
#outDir <- sprintf("%s/2017_PanCancer_Survival/clinNets_170430/featSelNets", 
#		rootDir)
outDir <- sprintf("%s/2017_PanCancer_Survival/pathwaysOnly_170502",rootDir)
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

source("writeConsensusNets_oneSet.R")
dt <- format(Sys.Date(),"%y%m%d")

logFile <- sprintf("%s/getConsNets_pathways_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
	for (curSet in "KIRC") { # names(dirList)){
		cat("-----------------------------\n")
		cat(sprintf("%s\n",curSet))
		cat("-----------------------------\n")
		for (cutoff in 10) {
		cat("-----------------------------\n")
		cat(sprintf("\tcutoff = %i\n", cutoff))
		cat("-----------------------------\n")

	pctPass <- 0.7
		outPfx		<-sprintf("%s/%s",outDir,curSet)
		outPfx	<- writeConsensusNets(datDir=dirList[[curSet]],
			outPfx=outPfx,consCutoff=cutoff,pctPass=pctPass)

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
