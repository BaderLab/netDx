# test effect of single-vs-multiple datasets to judge performance
# improvement

require(netDx)
inFile <- "/home/spai/BaderLab/DREAM_AD/input/GuanLab/C3/ADNI_Training_Q3_new.csv_matlab.csv"

outRoot <- "/home/spai/BaderLab/DREAM_AD/output"
numCores <- 8

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
# data preparation

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/integrate_keepMMSE%s",outRoot,dt)
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

dat <- read.delim(inFile,sep=",",h=F,as.is=T)
colnames(dat)[1:7] <- c("ID","AGE","IS_MALE","PTEDUCAT",
	"APOE4","MMSE","STATUS")
colnames(dat)[8:ncol(dat)] <- paste("MRI_",8:ncol(dat),sep="")
rownames(dat) <- dat[,1]


# clinical/genetic networks have similarity by normalized distance.
alldat <- t(dat[,-c(1,7)])
# remove MMSE
#idx <- which(rownames(alldat)%in% "MMSE")
#alldat <- alldat[-idx,]

pheno <- dat[,c("ID","STATUS")]
tmp <- rep(NA,nrow(pheno))
tmp[which(pheno$STATUS==1)] <- "CN"
tmp[which(pheno$STATUS==2)] <- "LMCI"
tmp[which(pheno$STATUS==4)] <- "AD"
pheno$STATUS <- tmp # recode categories as such.

rm(dat)

subtypes <- unique(pheno$STATUS)
netSets <- list()
for (k in 1:nrow(alldat)) 
		netSets[[rownames(alldat)[k]]] <- rownames(alldat)[k]

### create PSN - one net per variable
cat("* Building networks\n")
netDir <- sprintf("%s/networks", outDir)
t0 <- Sys.time()
netList <- makePSN_NamedMatrix(alldat,rownames(alldat),netSets,netDir,
			simMetric="custom",customFunc=normDiff,sparsify=TRUE,
			verbose=TRUE,numCores=numCores)
print(Sys.time()-t0)

# now create database
cat("* Creating database\n")
t0 <- Sys.time()
dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)
print(Sys.time()-t0)

# ----------------------------------------------------------------
# data preparation

combList <- list(
	clinical=c("AGE","IS_MALE","PTEDUCAT","MMSE"),
	genetic=c("APOE4"),
	imaging=rownames(alldat)[grep("MRI_",rownames(alldat))],
	all="all"
)

for (rngSeed in seq(5,50,10)) {
	cat("-------------------------------------------\n")
	cat(sprintf("RNG seed = %i\n",rngSeed))
	cat("-------------------------------------------\n")
	curd <- sprintf("%s/rng%i",outDir,rngSeed)

	dir.create(curd)
	cat(sprintf("RNG %i\n", rngSeed))

	# split into train/test for performance classes
	pheno$TT_STATUS <- splitTestTrain(pheno,setSeed=rngSeed)
	curd <- sprintf("%s/rng%i",outDir,rngSeed)

	for (cur in names(combList)) {
		t0 <- Sys.time()
		outRes <- list()
		pDir <- sprintf("%s/%s",curd, cur)
		dir.create(pDir)

		for (g in subtypes) {
			qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
									 pheno$TT_STATUS%in% "TRAIN")]
			qFile <- sprintf("%s/%s_testQuery",pDir,g)
	
			# use only selected nets
			if (cur == "all") nets <- "all" 
			else nets <- paste(combList[[cur]],"_cont",sep="")
			GM_writeQueryFile(qSamps,nets,nrow(pheno),qFile)
	
			resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
			outRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),
										  pheno,g)
		}
		outClass <- GM_OneVAll_getClass(outRes)
		both <- merge(x=pheno,y=outClass,by="ID")
		save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
		write.table(both,file=sprintf("%s/predictionResults.txt",pDir),
					sep="\t",col=TRUE,row=FALSE,quote=FALSE)

		cat(sprintf("%s complete\n", cur))
		print(Sys.time()-t0)
	}
}
