# DREAM Alzheimer's Subchallenge 3 predictor.
# goal to predict diagnostic class by MR imaging
rm(list=ls())

inFile <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/clinical/ADNIMERGE.csv"
outRoot <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/output"
numCores <- 2

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
# prepare input data

dat <- read.delim(inFile,sep=",",h=T,as.is=T)
# keep only baseline
dat 	<- subset(dat, VISCODE %in% "bl") 	# baseline visit
dat <- dat[,c("PTID","AGE","DX_bl","APOE4","PTEDUCAT","Hippocampus_bl")]

# make status variable
colnames(dat)[1] <- "ID"
dat <- na.omit(dat)
dat <- dat[-which(dat$DX_bl == "SMC"),] # remove midlevel states
dat$DX_bl[which(dat$DX_bl %in% c("EMCI","LMCI"))] <- "MCI"
colnames(dat)[which(colnames(dat) %in% "DX_bl")] <- "STATUS"

pheno 	<- dat; 
dat		<- pheno[,-which(colnames(pheno) %in% c("ID","DX_bl","STATUS"))]; 
rownames(dat) <- pheno$ID
dat		<- t(dat)


# ----------------------------------------------------------------
# run classifier

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/Q3_%s",outRoot,dt)
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)
require(netDx)

logFile <- sprintf("%s/log.txt", outDir)
sink(logFile,split=TRUE)
tryCatch({

### create PSN
netDir <- sprintf("%s/networks", outDir)

# do numerical ones first
netSets <- list(
	"clin_educ"="PTEDUCAT",
	"clin_age_at_dx"="AGE",
	"gen_APOE4"="APOE4",
	"img_Hpc_volume"="Hippocampus_bl"
	)

netList <- makePSN_NamedMatrix(dat,rownames(dat),netSets,netDir,
			simMetric="custom",customFunc=normDiff,sparsify=FALSE,
			verbose=TRUE)

# create a single GM database with training and test samples
dbDir 	<- GM_createDB(netDir,pheno$ID,outDir,numCores=numCores)

subtypes <- unique(pheno$STATUS)

for (rngSeed in c(1)) {
	curd <- sprintf("%s/rng%i",outDir,rngSeed)
	dir.create(curd)
	cat(sprintf("RNG %i\n", rngSeed))

	# split into train/test for performance classes
	pheno$TT_STATUS <- splitTestTrain(pheno,setSeed=rngSeed)

	# datatype combinations to try.
	combList <- list(
		#genetic="APOE4",clinical=c("DX_bl","age_at_dx"),
		#imaging="Hpc_volume",
		all="all")

	for (cur in names(combList)) {
		t0 <- Sys.time()
		cat(sprintf("%s\n",cur)) 
		pDir <- sprintf("%s/%s",curd, cur)
		dir.create(pDir)
		outRes <- list()
		
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
		res <- merge(x=pheno,y=outClass,by="ID")
		save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
		write.table(res,file=sprintf("%s/predictionResults.txt",pDir),
					sep="\t",col=TRUE,row=FALSE,quote=FALSE)
		cat(sprintf("%s complete\n", cur))
		print(Sys.time()-t0)
	}
} # end rng loop

},error=function(ex){
	print(ex)
},finally={
	cat("Closing log\n")
	sink(NULL)
})






