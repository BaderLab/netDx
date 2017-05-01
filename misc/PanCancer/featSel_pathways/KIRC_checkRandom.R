#' KIRC - eval performance on random set of nets the same size as
#' consensus nets
rm(list=ls())
require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9


#rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
rootDir <- "/mnt/data2/BaderLab"
#inDir <- sprintf("%s/2017_TCGA_KIRC/input",rootDir)
#outRoot <- sprintf("%s/2017_TCGA_KIRC/output",rootDir)
inDir <- sprintf("%s/PanCancer_KIRC/input",rootDir)
outRoot <- sprintf("%s/PanCancer_KIRC/output",rootDir)

#### makeConsProfiles code block >>>>
consNetDir <- sprintf("%s/PanCancer_common/oneClinNet_featSelNets",
		rootDir)
netScoreFile <- list(
		SURVIVENO=sprintf("%s/KIRC__thresh10_pctPass1.00_SURVIVENO_netScores.txt",
			consNetDir),
		SURVIVEYES=sprintf("%s/KIRC__thresh10_pctPass1.00_SURVIVEYES_netScores.txt",
			consNetDir)
)

maxRng <- 25 				# num train/test splits to check performance on
consCutoff <- 10		# score a net must have to be a consensus net
consPctPass <- 0.7	# % of rounds it must have score >= consCutoff

analysisMode <- "randomNets" # none |consNets | bestConsNets | randomNets
# none = just generate nets of consensus profiles, don't run predictions
# consNets = predictions with GMdb of all nets
# bestConsNets = GM db with nets that have corr < 0.01
# randomNets = GMdb with randomly sampled nets of size consNets 

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

if (!analysisMode %in% c("none","consNets","bestConsNets","randomNets")){
	cat("invalid method")
}

# count the number of consensus nets per group and whether clinical
# is part of this. will affect sampling downstream.
netCount <- list()
incClin <- list()
if (analysisMode == "randomNets") {
for (nm in names(netScoreFile)) {
		netScores	<- read.delim(netScoreFile[[nm]],sep="\t",h=T,as.is=T)
		netNames 	<- netScores[,1]
		netScores <- netScores[,-1]
		tmp <- rowSums(netScores >= consCutoff)
		idx <- which(tmp >= floor(consPctPass*ncol(netScores)))
		cat(sprintf("%s: %i nets score >= %i for %i%% trials\n",
			nm, length(idx), consCutoff, round(consPctPass*100)))
		
		if ("clinical.profile" %in% netNames[idx]) {
			netCount[[nm]] <- length(idx)-1 # subtract 1 for clinical
			incClin[[nm]] <- TRUE				# force clin to be included
		} else 
			netCount[[nm]] <- length(idx)
			incClin[[nm]] <- FALSE			# when sampling random, don't include clin
}
}

for (sampRNG in seq(55,500,5)) {

if (!file.exists(outRoot)) dir.create(outRoot)
dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/%s_forceClin%s_%s",outRoot,analysisMode,sampRNG,dt)
#### <<<< makeConsProfiles code block 

# -----------------------------------------------------------
# process input
clinList <- list(clinical=c("age","grade","stage"))

inFiles <- list(
	clinical=sprintf("%s/KIRC_clinical_core.txt",inDir),
	survival=sprintf("%s/KIRC_binary_survival.txt",inDir),
	rna=sprintf("%s/KIRC_mRNA_core.txt",inDir)
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

dats <- list() #input data in different slots

# clinical
cat("\t* Clinical\n")
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1];
clinical$grade <- as.numeric(factor(clinical$grade))
clinical$stage <- as.numeric(factor(clinical$stage))
clinical$ID <- NULL
clinical <- t(clinical)
dats$clinical <- clinical; rm(clinical)

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
rm(pheno,pheno_nosurv)

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({
	outDir <- megaDir

	pheno_all$TT_STATUS <- "DUMMY" 
	
	# create nets 
	netDir <- sprintf("%s/networks",outDir) 
	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
	   path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)
	data(genes)
	
###	# run featsel once per subtype
	subtypes <- unique(pheno_all$STATUS)
	
	## ----class-prediction, eval=TRUE-------------------------
	# now create GM databases for each class
	# should contain train + test patients
	# and be limited to nets that pass feature selection
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		# get feature selected net names
		
#### makeConsProfiles code block >>>>
		## SP: change the filenames here
		pTally <- switch(analysisMode,
	       consNets={
			read.delim(
				sprintf("%s/KIRC__thresh%i_pctPass1.00_%s_consNets.txt", 
					consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
###		}, none={
###			read.delim(
###				sprintf("%s/KIRC_%s_consNets.txt", consNetDir,g),
###				sep="\t",h=T,as.is=T)[,1]
####		}, bestConsNets={
###			tmp <- read.delim(
###				sprintf("%s/KIRC_%s_consNets.txt", consNetDir,g),
###				sep="\t",h=T,as.is=T)[,1]
###			fName <- dir(consNetDir,pattern=sprintf("KIRC_%s_correlations",g))
###			bestSet <- read.delim(
###				sprintf("%s/%s",consNetDir,fName),sep="\t",h=T,as.is=T)
###			idx_tmp <- unlist(sapply(1:nrow(bestSet), 
###				function(k) max(bestSet[k,4:6])))
###			oldlen <- length(tmp)
###			tmp <- tmp[which(tmp %in% rownames(bestSet)[which(idx_tmp>2)])]
###			cat(sprintf("Filtered %i to %i best\n", oldlen,length(tmp)))
###			tmp
		}, randomNets= {
			pTally	<- read.delim(
				sprintf("%s/KIRC__thresh%i_pctPass1.00_%s_AllNets.txt", 
						consNetDir,consCutoff,g),sep="\t",h=T,as.is=T)[,1]
			tmp <- which(pTally == "clinical.profile")
			if (any(tmp)) pTally <- pTally[-tmp]
			cat(sprintf("**** Sampling %i nets randomly ****\n",netCount[[g]]))
			set.seed(sampRNG);
			pTally <- sample(pTally, netCount[[g]], replace=FALSE)
			# only add clinical by force if it was part of the consensus
			# nets for the group
			if (incClin[[g]]) pTally <- c("clinical.profile", pTally)
			cat(sprintf("Nets for this round:\n%s\n", 
				paste(pTally,collapse="\n")))

				pTally
		},{
			stop(sprintf("Invalid analysisMode = %s\n", analysisMode))
		})
#### <<< makeConsProfiles code block 
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		cat(sprintf("%s: %i pathways\n",g,length(pTally)))
		netDir <- sprintf("%s/networks",pDir)
		if (!file.exists(netDir)) dir.create(netDir,recursive=TRUE)
	
    #prepare nets for new db
    # RNA 
    idx <- which(names(pathwayList) %in% pTally)
    if (any(idx)) {
        cat(sprintf("RNA: included %i nets\n", length(idx)))
        tmp <- makePSN_NamedMatrix(dats$rna, rownames(dats$rna),
    			pathwayList[idx],netDir,verbose=F,
					numCores=numCores, writeProfiles=TRUE)
    }
            
   # clinical
		if ("clinical.profile" %in% pTally) {
       	netList2 <- makePSN_NamedMatrix(dats$clinical, 
				rownames(dats$clinical),clinList[idx],
           netDir,writeProfiles=TRUE,
           verbose=TRUE,numCores=numCores,append=TRUE)
		}
        
### makeConsProfiles code block >>>>
		if (!analysisMode %in% "none") {
		# create db
		dbDir <- GM_createDB(netDir,pheno$ID,pDir,numCores=numCores)

		cat("Running test queries\n")
		GMdir <- sprintf("%s/GM_results",pDir)
		dir.create(GMdir)
		for (rngNum in 1:maxRng) {
			cat(sprintf("(%i)", rngNum))	
			pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
										  setSeed=rngNum*5)
			# query of all training samples for this class
			qSamps <- pheno$ID[which(pheno_all$STATUS %in% g & 
								 pheno_all$TT_STATUS%in%"TRAIN")]
	
			qFile <- sprintf("%s/%s_query_%i",GMdir,g,rngNum)
			GM_writeQueryFile(qSamps,"all",nrow(pheno_all),qFile)
			resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=GMdir)
			if (length(predRes)< rngNum) predRes[[rngNum]] <- list()
			predRes[[rngNum]][[g]] <- GM_getQueryROC(
									sprintf("%s.PRANK",resFile),
									pheno_all,g)
			unlink(resFile)
		}
		}
	}

	resDir <- sprintf("%s/predictions",outDir)
	dir.create(resDir)
	for (k in 1:length(predRes)) {
		predClass <- GM_OneVAll_getClass(predRes[[k]])
 		out <- merge(x=pheno_all,y=predClass,by="ID")
		outFile <- sprintf("%s/predictionResults_%i.txt",resDir,k)
	 	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
	}
### <<<< makeConsProfiles code block 
}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})
}
