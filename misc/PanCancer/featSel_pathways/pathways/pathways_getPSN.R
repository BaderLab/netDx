#' generate integrated PSN for TCGA breast cancer data
rm(list=ls())
# --------------------------------------------------------------
# data dirs for input
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
outDir <- sprintf("%s/2017_PanCancer_Survival/pathwaysOnly_170502",rootDir)
dataDir <- sprintf("%s/2017_TCGA_KIRC/output/pathOnly_consNets_170509",
	rootDir)
inData <- list(
setName="pathways",
clinical_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",
	rootDir),
survival_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",
	rootDir),
outDir=outDir,
ptFile=list(
	YES=sprintf("%s/SURVIVEYES/tmp/GENES.txt",dataDir),
	NO=sprintf("%s/SURVIVENO/tmp/GENES.txt",dataDir)
),
# net ID-to-name mappings
netInfo=list(
	YES=sprintf("%s/SURVIVEYES/tmp/NETWORKS.txt",dataDir),
	NO=sprintf("%s/SURVIVENO/tmp/NETWORKS.txt",dataDir)
),
# interaction nets
netDir=list(
	YES=sprintf("%s/SURVIVEYES/tmp/INTERACTIONS",dataDir),
	NO=sprintf("%s/SURVIVENO/tmp/INTERACTIONS",dataDir)
),
# we are going to take union of FS pathways for each class so we need
# the pathway scores for each class
netScoreFile=list(
	YES=sprintf("%s/featSelNets/KIRC_thresh10_pctPass0.70_SURVIVEYES_netScores.txt",
	outDir),
	NO=sprintf("%s/featSelNets/KIRC_thresh10_pctPass0.70_SURVIVENO_netScores.txt",
	outDir)
)
)
# --------------------------------------------------------------
# Work begins
logFile <- sprintf("%s/getPSN.log",outDir)
sink(logFile,split=TRUE)

tryCatch({
		source("../../getPSN_doItAll.R")
		getPSN(infoList=inData,consCutoff=10,consPctPass=0.7,topX=0.2,
			outDir=outDir)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		

