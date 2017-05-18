#' generate consensus maps for all PanCancer feature selection
#' predictors
require(netDx)
require(netDx.examples)

# root directory where all pancancer information is stored on the running
# machine
dirBase <- "/Users/shraddhapai/Documents/Research/BaderLab"
outDir <- sprintf("%s/2017_PanCancer_Survival/consNetsTest", dirBase)
dt <- format(Sys.Date(),"%y%m%d")

# with feature selection
dirList <- list(
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/output/OV_170227",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",
		dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222",dirBase)
	)

# directories with input files
inDirs <- list( 
	GBM=sprintf("%s/2017_TCGA_GBM/input",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/input",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/input",dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/input",dirBase)
)

# common data sources
# RNA-based pathway nets
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

#### ------------------------------------------------------------------------
### KIRC
### load genes in this dataset, to create an input set that only has genes
### measured here.

logFile <- sprintf("%s/KIRC_callConsMapMaker_%s.log",outDir,dt)
sink(logFile,split=TRUE)
tryCatch({
# somatic mutations
mutFile		<- sprintf("%s/KIRC_core_somatic_mutations.txt",inDirs$KIRC)
mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]
mutList <- pathwayList; 
names(mutList) <- paste("MUT_",names(mutList),sep="")
mutList <- lapply(mutList, function(x) x[which(x %in% mut_genes)])
ln <- unlist(lapply(mutList, length))
if (any(ln<1)) mutList <- mutList[-which(ln<1)]

# gene expression
xprFile <- sprintf("%s/KIRC_mRNA_core.txt",inDirs$KIRC)
xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
bpos <- regexpr("\\|", xpr_genes)
xpr_genes <- substr(xpr_genes, 1,bpos-1)
xprList <- pathwayList
xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])

netInfo <- list(
	mutations=list(pattern="MUT_", netSet=mutList),
	clinical=list(pattern="age|grade|stage",netSet="base"),
	RNA=list(pattern="*", netSet=xprList)
	)

cat("got here\n")

source("PanCancer_featSel_writeConsensusMap.R")
writeConsensusMap(datDir=dirList$KIRC, 
	predClasses=c("SURVIVEYES","SURVIVENO"),consCutoff=6L,
	netInfo=netInfo,outPfx=sprintf("%s/KIRC_%s_", outDir,dt))
rm(netInfo)
}, error=function(ex) {
	print(ex)
},finally={
	sink(NULL)
})

# ------------------------------------------------------------------------
###### GBM
###logFile <- sprintf("%s/GBM_callConsMapMaker_%s.log",outDir,dt)
###sink(logFile,split=TRUE)
###tryCatch({
#### somatic mutations
###mutFile		<- sprintf("%s/GBM_core_somatic_mutations.txt",inDirs$GBM)
###mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]
###mutList <- pathwayList; 
###names(mutList) <- paste("MUT_",names(mutList),sep="")
###mutList <- lapply(mutList, function(x) x[which(x %in% mut_genes)])
###ln <- unlist(lapply(mutList, length))
###if (any(ln<1)) mutList <- mutList[-which(ln<1)]
###
#### gene expression
###xprFile <- sprintf("%s/GBM_mRNA_core.txt",inDirs$GBM)
###xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
###xpr_genes <- sub("mRNA_","",xpr_genes)
###xprList <- pathwayList
###xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])
###
###netInfo <- list(
###	mutations=list(pattern="MUT_", netSet=mutList),
###	clinical=list(pattern="age|gender|Karnofsky",netSet="base"),
###	RNA=list(pattern="*", netSet=xprList)
###	)
###
###source("PanCancer_featSel_writeConsensusMap.R")
###writeConsensusMap(datDir=dirList$GBM, 
###	predClasses=c("SURVIVEYES","SURVIVENO"),
###	netInfo=netInfo,outPfx=sprintf("%s/GBM_%s_",outDir,dt),
###	consCutoff=6L)
###rm(netInfo)
###},error=function(ex) {
###	print(ex)
###},finally={
###	sink(NULL)
###})

# ------------------------------------------------------------------------
###### OV
###### load genes in this dataset, to create an input set that only has genes
###### measured here.
###
###logFile <- sprintf("%s/OV_callConsMapMaker_%s.log",outDir,dt)
###sink(logFile,split=TRUE)
###tryCatch({
#### somatic mutations
###mutFile		<- sprintf("%s/OV_core_somatic_mutations.txt",inDirs$OV)
###mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]
###mutList <- pathwayList; 
###names(mutList) <- paste("MUT_",names(mutList),sep="")
###mutList <- lapply(mutList, function(x) x[which(x %in% mut_genes)])
###ln <- unlist(lapply(mutList, length))
###if (any(ln<1)) mutList <- mutList[-which(ln<1)]
###
#### gene expression
###xprFile <- sprintf("%s/OV_mRNA_core.txt",inDirs$OV)
###xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
###xpr_genes <- sub("mRNA_","",xpr_genes)
###xprList <- pathwayList
###xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])
###
###netInfo <- list(
###	mutations=list(pattern="MUT_", netSet=mutList),
###	clinical=list(pattern="age|stage",netSet="base"),
###	protein=list(pattern="NOTCH3", netSet="base"),
###	DNAm=list(pattern="BRCA[12]", netSet="base"),
###	RNA=list(pattern="*", netSet=xprList)
###	)
###
###cat("got here\n")
###
###source("PanCancer_featSel_writeConsensusMap.R")
###writeConsensusMap(datDir=dirList$OV, 
###									predClasses=c("SURVIVEYES","SURVIVENO"),
###									netInfo=netInfo,
###									outPfx=sprintf("%s/OV_%s_", outDir,dt))
###rm(netInfo)
###}, error=function(ex) {
###	print(ex)
###},finally={
###				sink(NULL)
###})

#### ------------------------------------------------------------------------
###### LUSC 
###### load genes in this dataset, to create an input set that only has genes
###### measured here.
###
###logFile <- sprintf("%s/LUSC_callConsMapMaker_%s.log",outDir,dt)
###sink(logFile,split=TRUE)
###tryCatch({
#### somatic mutations
###mutFile		<- sprintf("%s/LUSC_core_somatic_mutations.txt",inDirs$LUSC)
###mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]
###mutList <- pathwayList; 
###names(mutList) <- paste("MUT_",names(mutList),sep="")
###mutList <- lapply(mutList, function(x) x[which(x %in% mut_genes)])
###ln <- unlist(lapply(mutList, length))
###if (any(ln<1)) mutList <- mutList[-which(ln<1)]
###
#### gene expression
###xprFile <- sprintf("%s/LUSC_mRNA_core.txt",inDirs$LUSC)
###xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
###xpr_genes <- sub("mRNA_","",xpr_genes)
###xprList <- pathwayList
###xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])
###
###netInfo <- list(
###	mutations=list(pattern="MUT_", netSet=mutList),
###	protein=list(pattern="PROT_", netSet="base"),
###	clinical=list(pattern="age|stage",netSet="base"),
###	RNA=list(pattern="*", netSet=xprList)
###	)
###
###source("PanCancer_featSel_writeConsensusMap.R")
###writeConsensusMap(datDir=dirList$LUSC, 
###									predClasses=c("SURVIVEYES","SURVIVENO"),
###									netInfo=netInfo,
###									outPfx=sprintf("%s/LUSC_%s_", outDir,dt))
###rm(netInfo)
###}, error=function(ex) {
###	print(ex)
###},finally={
###				sink(NULL)
###})
