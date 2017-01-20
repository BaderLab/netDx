# test multidata integration using ADNIMERGE data
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

# compute delta MMSEboth <- merge(x=bl,y=m24,by="PTID")
bl 	<- subset(dat, VISCODE %in% "bl") 	# baseline visit
m24	<- subset(dat, VISCODE %in% "m24")	# 24 month visit
both<- merge(x=bl,y=m24,by="PTID")
both$deltaMMSE <- both$MMSE.y-both$MMSE_bl.x

# keep only variables we will use
both <- both[,c("PTID","DX_bl.x","AGE.x","MMSE_bl.x","MMSE.y",
				"deltaMMSE","APOE4.x","Hippocampus_bl.x")]

# remove those with missing phenotype
idx <- which(is.na(both$deltaMMSE))
if (any(idx)) both <- both[-idx,]

# TODO remove those with delta_MMSE=-1 because these may not be 
# significantly worse.

# make status variable
both$DX <- factor(both$DX_bl.x, levels=c("CN","SMC","EMCI","LMCI","AD"))
both$STATUS <- "same"
both$STATUS[which(both$deltaMMSE< -1)] <- "worse"
colnames(both)[1] <- "ID"
both <- na.omit(both)
both <- both[-which(both$DX_bl.x == "SMC"),] # remove midlevel states
both$DX_bl.x[which(both$DX_bl.x %in% c("EMCI","LMCI"))] <- "MCI"

pheno 	<- both; 
dat		<- pheno[,-which(colnames(pheno) %in% c("ID","DX_bl.x","STATUS"))]; 
rownames(dat) <- pheno$ID
dat$DX <- as.integer(factor(pheno$DX_bl.x, levels=c("CN","MCI","AD")))
dat		<- t(dat)

rm(both,bl,m24,idx)

#### test similarity
###blah <- both$AGE.x
###names(blah) <- both$PTID
###sim <- normDiff(blah)
###require(reshape2)
###sim2 <- melt(sim)
###sim2 <- na.omit(sim2)
###
###aged <- dist(blah)
###aged <- melt(aged)

# ----------------------------------------------------------------
# run classifier

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/run%s_noSMC_MCI",outRoot,dt)
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
	"age_at_dx"="AGE.x",
	"APOE4"="APOE4.x",
	"Hpc_volume"="Hippocampus_bl.x",
	"DX_bl"="DX"
	)

netList <- makePSN_NamedMatrix(dat,rownames(dat),netSets,netDir,
			simMetric="custom",customFunc=normDiff,sparsify=FALSE,
			verbose=TRUE)

# create a single GM database with training and test samples
dbDir 	<- GM_createDB(netDir,pheno$ID,outDir,numCores=numCores)

subtypes <- unique(pheno$STATUS)

for (rngSeed in c(1:10)) {
	curd <- sprintf("%s/rng%i",outDir,rngSeed)
	dir.create(curd)
	cat(sprintf("RNG %i\n", rngSeed))

	# split into train/test for performance classes
	pheno$TT_STATUS <- splitTestTrain(pheno,setSeed=rngSeed)

	# datatype combinations to try.
	combList <- list(
		genetic="APOE4",clinical=c("DX_bl","age_at_dx"),
		imaging="Hpc_volume",
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
		both <- merge(x=pheno,y=outClass,by="ID")
		save(outRes,file=sprintf("%s/outRes.Rdata",pDir))
		write.table(both,file=sprintf("%s/predictionResults.txt",pDir),
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






