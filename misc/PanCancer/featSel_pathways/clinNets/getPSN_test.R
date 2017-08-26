#' generate integrated PSN for TCGA breast cancer data
rm(list=ls())

# --------------------------------------------------------------
# Param for computing integrated PSN
consCutoff 		<-10  	# include nets with score >= this value
consPctPass		<- 1
dt <- format(Sys.Date(),"%y%m%d")
topX <- 0.2
aggFun <- "MEAN"
setName <- "clinNets"

# --------------------------------------------------------------
# data dirs for input

#### CHANGE THIS TO THE ROOTDIR WHERE YOU UNZIP KIRC_PSNtest
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab/KIRC_PSNtest"
####

dataDir	<- sprintf("%s/output/KIRC_clinNets_170430", rootDir)
outDir 	<- sprintf("%s/results",rootDir)
inData <- list(
setName="clinNets",
clinical_file=sprintf("%s/input/KIRC_clinical_core.txt",
	rootDir),
survival_file=sprintf("%s/input/KIRC_binary_survival.txt",
	rootDir),
outDir=outDir,

### Points to GeneMANIA db within one of the rngs/ to get a list 
### of patients, network names, and the interaction nets (which have been
### renamed to internal ids, such as 1.1.txt.

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
	YES=sprintf("%s/KIRC__thresh10_pctPass1.00_SURVIVEYES_netScores.txt",
outDir),
	NO=sprintf("%s/KIRC__thresh10_pctPass1.00_SURVIVENO_netScores.txt",
outDir)
)
)

# --------------------------------------------------------------
# Work begins
logFile <- sprintf("%s/getPSN.log",outDir)
sink(logFile,split=TRUE)

tryCatch({
		source("../../getPSN_doItAll.R")

		### note that right now the function is taking a lot of input from
		### the predictor. This is because in a complex design where different
		### nets may be created in different ways, it may be easier to have
		### another function create the nets and just pass it to this function.
		getPSN(infoList=inData,consCutoff=10,consPctPass=1,topX=0.2,
			outDir=outDir)
},error=function(ex){ 
	print(ex)
}, finally={
	sink(NULL)
})
		
