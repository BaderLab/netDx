# --------------------------------------------------------------
# writes high-scoring nets for each cancer type
rm(list=ls())
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_pathway_clinOneNet_170428",
		rootDir)
)

source("writeConsensusNets_oneSet.R")
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/2017_PanCancer_Survival/oneClinNet_featSelNets", 
		rootDir)
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

pctPass <- 1.0 

logFile <- sprintf("%s/getConsNets_pathways_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
	for (curSet in "KIRC") { # names(dirList)){
		cat("-----------------------------\n")
		cat(sprintf("%s\n",curSet))
		cat("-----------------------------\n")
		for (cutoff in 9:10) {
		cat("-----------------------------\n")
		cat(sprintf("\tcutoff = %i; pct pass = %1.2f\n", cutoff,pctPass))
		cat("-----------------------------\n")

		outPfx<-sprintf("%s/%s_",outDir,curSet,cutoff)
		writeConsensusNets(datDir=dirList[[curSet]],
			outPfx=outPfx,consCutoff=cutoff,pctPass=pctPass)

		cat("YES\n")
		fname <- sprintf("%s_thresh%i_pctPass%1.2f_SURVIVEYES_consNets.txt",
			outPfx,cutoff,pctPass)
		dat <- read.delim(fname,sep="\t",header=T,as.is=T)
		cat("\n")
		cat("NO\n")
		fname <- sprintf("%s_thresh%i_pctPass%1.2f_SURVIVENO_consNets.txt",
			outPfx,cutoff,pctPass)
		dat <- read.delim(fname,sep="\t",header=T,as.is=T)
		print(dat)

		cat("-----------------------------\n")
		}
	}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
