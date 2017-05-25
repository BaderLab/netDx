#' generate integrated PSN for KIRC - oneNetPer.
rm(list=ls())
# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
dataDir <- sprintf("%s/2017_TCGA_KIRC/output/KIRC_oneNetPer_normDiff_170518", 
	rootDir)
outDir <- sprintf("%s/2017_PanCancer_Survival/oneNetPer_FeatSel",rootDir)
inData <- list(
	setName="oneClinNet",
clinical_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",
	rootDir),
survival_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",
	rootDir),
outDir=outDir,
ptFile=list(YES=sprintf("%s/tmp/GENES.txt",dataDir)),
# net ID-to-name mappings
netInfo=list(YES=sprintf("%s/tmp/NETWORKS.txt",dataDir)),
# interaction nets
netDir=list(YES=sprintf("%s/tmp/INTERACTIONS",dataDir))
)
# --------------------------------------------------------------
# Work begins

logFile <- sprintf("%s/getPSN.log",outDir)
sink(logFile,split=TRUE)

tryCatch({
		source("../getPSN_doItAll.R")
		getPSN(infoList=inData,consCutoff=10,consPctPass=1,topX=0.5)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
