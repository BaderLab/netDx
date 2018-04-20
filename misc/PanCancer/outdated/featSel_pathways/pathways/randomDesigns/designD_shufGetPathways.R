#' KIRC random design D get pathways.
rm(list=ls())
require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9

rootDir <- "/home/shraddhapai/BaderLab"
inDir <- sprintf("%s/PanCancer_KIRC/input",rootDir)
outRoot <- sprintf("%s/PanCancer_KIRC/output",rootDir)

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/pathRandom_designD_Dummy_%s",outRoot,dt)

# -----------------------------------------------------------
## count the number of consensus nets per group and whether clinical
### is part of this. will affect sampling downstream.
consNetDir <- sprintf("%s/PanCancer_common/pathwaysOnly_170502",rootDir)
netScoreFile <- list(
	SURVIVENO=sprintf("%s/pathways_thresh10_pctPass0.70_SURVIVENO_netScores.txt",
		consNetDir),
	SURVIVEYES=sprintf("%s/pathways_thresh10_pctPass0.70_SURVIVEYES_netScores.txt",
		consNetDir)
)
netCount <- list()
for (nm in names(netScoreFile)) {
   netScores	<- read.delim(netScoreFile[[nm]],sep="\t",h=T,as.is=T)
   netNames 	<- netScores[,1]
   netScores <- netScores[,-1]
   netCount[[nm]] <- colSums(netScores>=cutoff,na.rm=TRUE)
   cat(sprintf("%s: # fs nets\n", nm))
   print(summary(netCount[[nm]]))
}

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"
#======transform clinical data=========
pheno$grade <- as.vector(pheno$grade)
pheno$grade[pheno$grade=="G1"] <- "G2"
pheno$grade[pheno$grade=="GX"] <- "G2"
pheno$grade <- as.factor(pheno$grade)
pheno <- pheno[, -which(colnames(pheno)=="gender")]
#======================================

surv <- read.delim(inFiles$survival,sep="\t",h=T,as.is=T)
colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")
pheno$X <- NULL
# pheno$gender <- ifelse(pheno$gender=="FEMALE",1, 0)
pheno_nosurv <- pheno[1:4]

pheno_all <- pheno; 
rm(pheno,pheno_nosurv)
# ----------------------------------------------------------
# build classifier
dir.create(megaDir)

# list networks
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
netFile <- sprintf("%s/inputNets.txt", megaDir)

tryCatch({
for (rngNum in 1:100) {
	# to force rng to set
	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
											  setSeed=rngNum*5)
	pheno <- pheno_all
	subtypes <- unique(pheno_all$STATUS)
	predRes <- list()
	for (g in subtypes) {

		# randomly sample net names
		pTally	<- read.delim(
			sprintf("%s/pathways_thresh10_pctPass0.70_%s_AllNets.txt", 
				consNetDir,g),sep="\t",h=T,as.is=T)[,1]
		pTally <- sub(".profile","",pTally)
		cur_randNum <- netCount[[g]][rngNum]
		pTally <- sample(pTally, cur_randNum, replace=FALSE) ## random
			
		write.table(pTally,
			file=sprintf("%s/nets_%s_%i.txt",megaDir,g,rngNum),
			sep="\t",col=F,row=F,quote=F)
	}
}
}, error=function(ex){
	print(ex)
}, finally={
})

