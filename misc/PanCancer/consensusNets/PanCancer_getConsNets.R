# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222",rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",
		rootDir),
	OV=sprintf("%s/2017_TCGA_OV/output/OV_170227",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",
		rootDir)
)

source("PanCancer_featSel_writeConsensusNets.R")
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/2017_PanCancer_Survival/consNetsTest", rootDir)
if (!file.exists(outDir)) dir.create(outDir)
logFile <- sprintf("%s/getConsNets_%s.log",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
for (curSet in names(dirList)){
	cat("-----------------------------\n")
	cat(sprintf("%s\n",curSet))
	cat("-----------------------------\n")
	for (cutoff in 6:9) {
	cat(sprintf("\tcutoff = %i\n", cutoff))
	writeConsensusNets(datDir=dirList[[curSet]],
		outPfx=sprintf("%s/%s_thresh%i",outDir,curSet,cutoff),
		consCutoff=cutoff)
	}
}
},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
