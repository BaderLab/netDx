#' feature selection for GBM from PanCancer survival dataset
#' 10-fold CV predictor design with clinical and mRNA data
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
megaDir <- sprintf("%s/pathways_AltClass_%s",outRoot,dt)

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

netScoreList <- list()
for (nm in names(netScoreFile)) {
	netScoreList[[nm]] <- list()
  netScores	<- read.delim(netScoreFile[[nm]],sep="\t",h=T,as.is=T)
  netNames 	<- netScores[,1]
  netScores <- netScores[,-1]
	for (k in 1:ncol(netScores)) { 
		netScoreList[[nm]][[k]] <- netNames[which(netScores[,k]>=9)]
		browser()
	}
}

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir),
	rna=sprintf("%s/KIRC_mRNA_core.txt",inDir),
	mut=sprintf("%s/from_firehose/KIRC_core_somatic_mutations.txt",
		inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

#======transform clinical data=========
# KIRC-specific
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

dats <- list() #input data in different slots

# RNA
cat("\t* RNA\n")
rna <- read.delim(inFiles$rna,sep="\t",h=T,as.is=T)
rna <- t(rna)
colnames(rna) <- rna[1,]; rna <- rna[-1,];
rna <- rna[-nrow(rna),]
class(rna) <- "numeric"
rownames(rna) <- sub("mRNA_","",rownames(rna))
rownames(rna) <- sub("\\..*","",rownames(rna))
dats$rna <- rna; rm(rna)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) {
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})

pheno_all <- pheno;
subtypes <- unique(pheno$STATUS)
rm(pheno,pheno_nosurv)

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

# list networks
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
netFile <- sprintf("%s/inputNets.txt", megaDir)
cat(sprintf("NetType\tNetName\n",file=netFile))
# rna nets
for (i in names(pathwayList)) {
	cat(sprintf("rna\t%s\n",i),file=netFile,append=TRUE)
}

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({
for (rngNum in 1:25) {
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))
	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)
	
	# no feature selection

	## ----class-prediction, eval=TRUE-------------------------
	# now create GM databases for each class
	# should contain train + test patients
	# and be limited to nets that pass feature selection
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		# get feature selected net names
		pTally <- netScoreList[[g]][[rngNum]]
		cat(sprintf("%s: rng %i: got %i nets\n", g, rngNum, length(pTally)))
		browser()
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		cat(sprintf("%s: %i pathways\n",g,length(pTally)))
		netDir <- sprintf("%s/networks",pDir)

     # prepare nets for new db
     # RNA
     idx <- which(names(pathwayList) %in% pTally)
     if (any(idx)) {
       cat(sprintf("RNA: included %i nets\n", length(idx)))
       tmp <- makePSN_NamedMatrix(dats$rna, rownames(dats$rna),
          pathwayList[idx],
         netDir,verbose=F,numCores=numCores, writeProfiles=TRUE)
    }

		### Alternate way of classifying patients
		### for each group, build a database only containing TRAIN for the 
		### class of interest and all TEST.
		### leave out TRAIN from the other classes.
		curPheno <- pheno
		idx1 <- which(pheno$TT_STATUS %in% "TEST")
		idx2 <- which(pheno$TT_STATUS %in% "TRAIN" & pheno$STATUS %in% g)
		cat(sprintf("Making db with %i %s TRAIN and %i TEST\n",
			length(idx2),g,length(idx1)))
		curPheno <- pheno[union(idx1,idx2),]

		# create db
		dbDir <- GM_createDB(netDir,curPheno$ID,pDir,numCores=numCores)
		# query of all training samples for this class
		qSamps <- curPheno$ID[which(curPheno$STATUS %in% g &
								 curPheno$TT_STATUS%in%"TRAIN")]

		qFile <- sprintf("%s/%s_query",pDir,g)
		GM_writeQueryFile(qSamps,"all",nrow(curPheno),qFile)
		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
		predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),curPheno,g)
	}

	predClass <- GM_OneVAll_getClass(predRes)
	out <- merge(x=pheno_all,y=predClass,by="ID")
	outFile <- sprintf("%s/predictionResults.txt",outDir)
	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)

	acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
	cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
		nrow(out), acc*100))

	require(ROCR)
	ROCR_pred <- prediction(out$SURVIVEYES_SCORE-out$SURVIVENO,
						out$STATUS=="SURVIVEYES")
	save(predRes,ROCR_pred,file=sprintf("%s/predRes.Rdata",outDir))

    #cleanup to save disk space
    system(sprintf("rm -r %s/dataset %s/tmp %s/networks",
        outDir,outDir,outDir))
    system(sprintf("rm -r %s/SURVIVENO/dataset %s/SURVIVENO/networks",
        outDir,outDir))
    system(sprintf("rm -r %s/SURVIVEYES/dataset %s/SURVIVEYES/networks",
        outDir,outDir))
}
}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})


