#' generate integrated PSN for TCGA breast cancer data
rm(list=ls())
require(netDx)
# --------------------------------------------------------------
# Param for computing integrated PSN
consCutoff 		<-10  	# include nets with score >= this value
consPctPass		<- 1
dt <- format(Sys.Date(),"%y%m%d")

netMode <- "consensus" # consensus|bestAUC
# --------------------------------------------------------------
# data dirs for input
rootDir	<- "/Users/shraddhapai/DropBox/netDx/BaderLab"
outDir	<- sprintf("%s/2017_PanCancer_Survival/clinRNA_best",
		rootDir)
dataDir <- sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinRNA_best", 
			rootDir)
setName <- "clinRNA_best"

inData 	<- list(
setName="clinRNA_best",
clinical_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",
	rootDir),
survival_file=sprintf("%s/2017_TCGA_KIRC/input/KIRC_binary_survival.txt",
	rootDir),
outDir=outDir,
ptFile=list(
	YES=sprintf("%s/rng1/SURVIVEYES/tmp/GENES.txt",dataDir),
	NO=sprintf("%s/rng1/SURVIVENO/tmp/GENES.txt",dataDir)
),
# net ID-to-name mappings
netInfo=list(
	YES=sprintf("%s/rng1/SURVIVEYES/tmp/NETWORKS.txt",dataDir),
	NO=sprintf("%s/rng1/SURVIVENO/tmp/NETWORKS.txt",dataDir)
),
# interaction nets
netDir=list(
	YES=sprintf("%s/rng1/SURVIVEYES/tmp/INTERACTIONS",dataDir),
	NO=sprintf("%s/rng1/SURVIVENO/tmp/INTERACTIONS",dataDir)
),
# we are going to take union of FS pathways for each class so we need
# the pathway scores for each class
netScoreFile=list(
	YES=sprintf("%s/clinRNA_best_thresh10_pctPass0.70_SURVIVEYES_netScores.txt",
			outDir),
	NO=sprintf("%s/clinRNA_best_thresh10_pctPass0.70_SURVIVENO_netScores.txt",
			outDir)
))

# --------------------------------------------------------------
# Work begins

logFile <- sprintf("%s/getPSN_%s.log",outDir,dt)
sink(logFile,split=TRUE)

tryCatch({
		source("../../getPSN_doItAll.R")
		getPSN(infoList=inData,consCutoff=10,consPctPass=1,topX=0.2,
			outDir=outDir)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		
