# --------------------------------------------------------------
# writes high-scoring nets for each cancer type
rm(list=ls())
source("writeConsensusNets_oneSet.R")

set2run <- "pathway80"
cutoff <- 10
pctPass <- 0.7
rnaFile <- "/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/input/KIRC_mRNA_core.txt"

# data dirs for input
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
dirList <- list(
	clinNets=sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinNets_170430",rootDir),
	clinRNA_best=sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinRNA_best",rootDir),
	pathways=sprintf("%s/2017_TCGA_KIRC/output/pathway_170502",rootDir),
	pathway80=sprintf("%s/2017_TCGA_KIRC/output/pathway80_170719",rootDir)
)

outList <- list(
		clinNets=sprintf("%s/2017_PanCancer_Survival/clinNets_170430/featSelNets",
		rootDir),
		pathways=sprintf("%s/2017_PanCancer_Survival/pathwaysOnly_170502/featSelNets",rootDir),
		pathway80=sprintf("%s/2017_PanCancer_Survival/pathway80_170502/featSelNets",rootDir)
)

outDir <- outList[[set2run]]
if (!file.exists(outDir)) dir.create(outDir,recursive=TRUE)

dt <- format(Sys.Date(),"%y%m%d")
logFile <- sprintf("%s/getConsNets_pathways_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
	rna <- read.delim(rnaFile,sep="\t",h=T,as.is=T)
	for (curSet in set2run){
		cat("-----------------------------\n")
		cat(sprintf("%s\n",curSet))
		cat("-----------------------------\n")

		outPfx		<-sprintf("%s/%s",outDir,curSet)
		outPfx	<- writeConsensusNets(datDir=dirList[[curSet]],
			outPfx=outPfx,consCutoff=10,pctPass=pctPass)

		cat("-----------------------------\n")
	}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
