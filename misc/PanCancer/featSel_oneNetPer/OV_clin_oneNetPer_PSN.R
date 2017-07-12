#' generate integrated PSN for OV - oneNetPer.
rm(list=ls())

# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
dataDir <- sprintf("%s/2017_TCGA_OV/output/OV_oneNetPer_170425", 
	rootDir)

outDir <- sprintf("%s/2017_PanCancer_Survival/oneNetPer_FeatSel/OV",rootDir)
inData <- list(
	setName="OV_oneClinNet",
	clinical_file=sprintf("%s/2017_TCGA_OV/input/OV_clinical_core.txt",
		rootDir),
	survival_file=sprintf("%s/2017_TCGA_OV/input/OV_binary_survival.txt",
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
		getPSN(infoList=inData,consCutoff=10,consPctPass=1,topX=0.7,outDir=outDir)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})

