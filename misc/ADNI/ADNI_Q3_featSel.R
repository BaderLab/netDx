#' feature selection for DREAM AD Question 3
#' here we start with the input files processed by GuanLab, one of the
#' winning teams of DREAM 2014 AD challenge.
#' 
#' we perform feature selection over the three diagnostic groups and 
#' then do a final classification.

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.67
cutoff <- 9

require(netDx)
inFile <- "/home/spai/BaderLab/DREAM_AD/input/GuanLab/C3/ADNI_Training_Q3_new.csv_matlab.csv"
outRoot <- "/home/spai/BaderLab/DREAM_AD/output"

IMG_START <- 8

# ----------------------------------------------------------------
# helper functions

# normalized difference 
# x is vector of values, one per patient (e.g. ages)
normDiff <- function(x) {
    #if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x,na.rm=T)-min(x,na.rm=T)
    
    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(x[i]-x[j])/(max(x)-min(x)))
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}

# ----------------------------------------------------------------

dt <- format(Sys.Date(),"%y%m%d")
outRoot <- sprintf("%s/featSel20RNG_noMMSE_%s",outRoot,dt)
if (file.exists(outRoot)) unlink(outRoot,recursive=TRUE)
dir.create(outRoot)

logFile <- sprintf("%s/log.txt",outRoot)
sink(logFile,split=TRUE)
tryCatch({

dat <- read.delim(inFile,sep=",",h=F,as.is=T)
colnames(dat)[1:7] <- c("ID","AGE","IS_MALE","PTEDUCAT",
	"APOE4","MMSE","STATUS")
colnames(dat)[8:ncol(dat)] <- paste("MRI_",8:ncol(dat),sep="")
rownames(dat) <- dat[,1]


# clinical/genetic networks have similarity by normalized distance.
alldat <- t(dat[,-c(1,7)])
# remove MMSE
idx <- which(rownames(alldat)%in% "MMSE")
alldat <- alldat[-idx,]

pheno <- dat[,c("ID","STATUS")]
tmp <- rep(NA,nrow(pheno))
tmp[which(pheno$STATUS==1)] <- "CN"
tmp[which(pheno$STATUS==2)] <- "LMCI"
tmp[which(pheno$STATUS==4)] <- "AD"
pheno$STATUS <- tmp # recode categories as such.

pheno_mega <- pheno
dat_mega <- alldat
subtypes <- unique(pheno$STATUS)
netSets <- list()
for (k in 1:nrow(dat_mega)) 
		netSets[[rownames(dat_mega)[k]]] <- rownames(dat_mega)[k]

rm(pheno,dat,alldat) 

# ---------------------------------------------------------
# feature selection
for (setSeed in seq(70,100,5)) {
	cat("-------------------------------------------\n")
	cat(sprintf("RNG seed = %i\n",setSeed))
	cat("-------------------------------------------\n")

	pheno_all <- pheno_mega
	dat_all <- dat_mega

	outDir <- sprintf("%s/RNG%i",outRoot,setSeed)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
										  setSeed=setSeed)
	print(table(pheno_all[,c("STATUS","TT_STATUS")]))
	
	idx		<- which(pheno_all$TT_STATUS %in% "TRAIN")
	pheno	<- pheno_all[idx,]
	dat 	<- dat_all[,idx]
	
	### create PSN - one net per variable
	netDir <- sprintf("%s/networks", outDir)
	netList <- makePSN_NamedMatrix(dat,rownames(dat),netSets,netDir,
				simMetric="custom",customFunc=normDiff,sparsify=TRUE,
				verbose=TRUE,numCores=numCores)
	
	# now create database
	dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)
	
	# run 10-fold cv per subtype
	for (g in subtypes) {
	    pDir <- sprintf("%s/%s",outDir,g)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
		dir.create(pDir)
		cat(sprintf("\n******\nSubtype %s\n",g))
		pheno_subtype <- pheno
		
		## label patients not in the current class as a residual
		pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
		## sanity check
		print(table(pheno_subtype$STATUS,useNA="always"))
		resDir    <- sprintf("%s/GM_results",pDir)
		## query for feature selection comprises of training 
		## samples from the class of interest
		trainPred <- pheno$ID[which(pheno$STATUS %in% g)]
		
		# Cross validation
		GM_runCV_featureSet(trainPred, resDir, dbDir$dbDir, 
			nrow(pheno_subtype),verbose=T, numCores=numCores,
			GMmemory=GMmemory)
	
		# patient similarity ranks
		prank <- dir(path=resDir,pattern="PRANK$")
		# network ranks
		nrank <- dir(path=resDir,pattern="NRANK$")
		cat(sprintf("Got %i prank files\n",length(prank)))
			
	    # Compute network score
		pTally		<- GM_networkTally(paste(resDir,nrank,sep="/"))
		head(pTally)
		# write to file
		tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
		write.table(pTally,file=tallyFile,sep="\t",col=T,row=F,quote=F)
	}
	
	## ----class-prediction, eval=TRUE-------------------------
	# now create GM databases for each class
	# should contain train + test patients
	# and be limited to nets that pass feature selection
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		# get feature selected net names
		pTally <- read.delim(
			sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir,g),
			sep="\t",h=T,as.is=T)
		pTally <- pTally[which(pTally[,2]>=cutoff),1]
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		cat(sprintf("%s: %i pathways\n",g,length(pTally)))
		netDir <- sprintf("%s/networks",pDir)

		# prepare nets for new db
		tmp <- makePSN_NamedMatrix(dat_all,rownames(dat_all),
			netSets[which(names(netSets)%in% pTally)],
			simMetric="custom",customFunc=normDiff,sparsify=TRUE,
			netDir,verbose=F,numCores=numCores)
		# create db
		dbDir <- GM_createDB(netDir,pheno$ID,pDir,numCores=numCores)
		# query of all training samples for this class
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
								 pheno$TT_STATUS%in%"TRAIN")]

		qFile <- sprintf("%s/%s_query",pDir,g)
		GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
		predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
	}
	## ----label-predictions, eval=TRUE,cache=FALSE-------------------
	save(predRes,file=sprintf("%s/predRes.Rdata",outDir))
	predClass <- GM_OneVAll_getClass(predRes)
	out <- merge(x=pheno_all,y=predClass,by="ID")
	outFile <- sprintf("%s/predClass.txt",outDir)
	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
	
	acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
	cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
		nrow(out), acc*100))
	
	rm(pheno_all,dat_all) # cleanup.
	system(sprintf("rm -r %s/tmp %s/dataset %s/networks", outDir,
		outDir,outDir))
}

},error=function(ex){
	print(ex)
},finally={
	sink(NULL)
})
