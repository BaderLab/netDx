#' feature selection for LUSC  from PanCancer survival dataset
#' 10-fold CV predictor design with clinical and proteomic data

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9

inDir <- "/mnt/data2/BaderLab/PanCancer_LUSC/input"
outRoot <-"/mnt/data2/BaderLab/PanCancer_LUSC/output"

dt <- format(Sys.Date(),"%y%m%d")
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

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/LUSC_clinical_core.txt",inDir),
	survival=sprintf("%s/LUSC_binary_survival.txt",inDir),
	rna=sprintf("%s/LUSC_mRNA_core.txt",inDir),
	rppa=sprintf("%s/LUSC_RPPA_core.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"

surv <- read.delim(inFiles$survival,sep="\t",h=T,as.is=T)
colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")

dats <- list() #input data in different slots

# clinical
cat("\t* Clinical\n")
clin <- pheno
# ######
# This section copied from main.R of syn1895966 and adapted to the current
# code
clin$stage <- as.vector(clin$stage)
clin$stage[clin$stage=="Stage IA"| clin$stage=="Stage IB"] <- "I"
clin$stage[clin$stage=="Stage IIA"| clin$stage=="Stage IIB"| clin$stage=="Stage II"] <- "II"
clin$stage[clin$stage=="Stage IIIA"| clin$stage=="Stage IIIB"] <- "III"
clin$stage <- as.factor(clin$stage)
clin <- clin[, -which(colnames(clin)=="gender")]
# ######
rownames(clin) <- clin[,1]; 
clin <- t(clin[,c("age","stage")])
clin[1,] <- as.integer(clin[1,])
clin[2,] <- as.integer(as.factor(clin[2,]))
class(clin) <- "numeric"
dats$clin <- clin; rm(clin)

# Proteomics
cat("\t* RPPA\n")
rppa <- read.delim(inFiles$rppa,sep="\t",h=T,as.is=T)
rppa <- t(rppa)
colnames(rppa) <- rppa[1,]; rppa <- rppa[-1,]; 
rppa <- rppa[-nrow(rppa),]
class(rppa) <- "numeric"
nm <- sub("RPPA_","",rownames(rppa))
dpos <- regexpr("\\.",nm)
nm <- substr(nm,1,dpos-1)
rownames(rppa) <- nm

dats$rppa <- rppa; rm(rppa)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})

pheno_all <- pheno; rm(pheno)

# ----------------------------------------------------------
# build classifier
require(netDx)
require(netDx.examples)
outDir <- sprintf("%s/featSel_%s",outRoot,dt)
numCores <- 8L
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

setSeed <- 102
pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
										  setSeed=setSeed)

# --------------------------------------------
# feature selection - train only
pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
dats_train <- lapply(dats,function(x) { x[,which(colnames(x) %in% pheno$ID)]})

# create nets 
netDir <- sprintf("%s/networks",outDir) 
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
data(genes)
gene_GR     <- GRanges(genes$chrom,IRanges(genes$txStart,genes$txEnd),
  	name=genes$name2)
cat("* Limiting to pathway genes\n")
path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)

# group by pathway
netList <- makePSN_NamedMatrix(dats_train$rppa, rownames(dats_train$rppa),
        pathwayList,netDir,verbose=FALSE,                                      
        numCores=numCores,writeProfiles=TRUE)                                   
cat("Made protein pathway nets\n")

# each clinical var is its own net
clinList <- list(age="age",stage="stage")
netList2 <- makePSN_NamedMatrix(dats_train$clin, rownames(dats_train$clin),
		clinList,netDir, simMetric="custom",customFunc=normDiff,sparsify=TRUE,
		verbose=TRUE,numCores=numCores,append=TRUE)
cat("Made clinical nets\n")

netList <- unlist(c(netList,netList2)) 
cat(sprintf("Total of %i nets\n", length(netList)))

# now create database
dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)

# run featsel once per subtype
subtypes <- unique(pheno$STATUS)
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
	trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% g)]
	
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
	tmp <- makePSN_NamedMatrix(dats$rppa, rownames(dats$rppa),
        pathwayList[which(names(pathwayList)%in%pTally)],
		netDir,verbose=F,numCores=numCores)

	netList2 <- makePSN_NamedMatrix(dats$clin, rownames(dats$clin),
		clinList,netDir, simMetric="custom",customFunc=normDiff,
		sparsify=TRUE,
		verbose=TRUE,numCores=numCores,append=TRUE)

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

save(predRes,file=sprintf("%s/predRes.Rdata",outDir))
predClass <- GM_OneVAll_getClass(predRes)
out <- merge(x=pheno_all,y=predClass,by="ID")
outFile <- sprintf("%s/predClass.txt",outDir)
write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)

acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
	nrow(out), acc*100))

