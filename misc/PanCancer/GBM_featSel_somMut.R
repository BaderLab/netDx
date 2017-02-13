#' feature selection for GBM from PanCancer survival dataset
#' 10-fold CV predictor design with clinical and mRNA data
rm(list=ls())
require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9

inDir <- "/mnt/data2/BaderLab/PanCancer_GBM/input"
outRoot <-"/mnt/data2/BaderLab/PanCancer_GBM/output"
#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/input"
#outRoot <-"/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/output"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/featSel_%s",outRoot,dt)

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
	clinical=sprintf("%s/GBM_clinical_core.txt",inDir),
	survival=sprintf("%s/GBM_binary_survival.txt",inDir),
	rna=sprintf("%s/GBM_mRNA_core.txt",inDir),
	mut=sprintf("%s/from_firehose/GBM_core_somatic_mutations.txt",
		inDir)
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
pheno$X <- NULL
pheno$gender <- ifelse(pheno$gender=="FEMALE",1, 0)
pheno_nosurv <- pheno[1:4]

dats <- list() #input data in different slots

# clinical
cat("\t* Clinical\n")
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1]; 
clinical$ID <- NULL
clinical$performance_score[which(clinical$performance_score == "[Not Available]")] <- NA
clinical$performance_score <- strtoi(clinical$performance_score)
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
dats$rna <- rna; rm(rna)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})

# somatic mutations
cat("\t* Somatic mutations\n")
mut <- read.delim(inFiles$mut,sep="\t",h=T,as.is=T)
# next steps: convert into GRanges with LOCUS_NAMES column.
# call makePSN_NamedRanges() to make pathway-level nets.
# run training.
pat_GR <- GRanges(paste("chr",mut$Chromosome,sep=""),
	IRanges(mut$Start_position,mut$End_position))
pat_GR$LOCUS_NAMES<-mut$Hugo_Symbol
pat_GR$ID <- mut$ID

pheno_all <- pheno; 
pat_GR_all <- pat_GR;

rm(pheno,pheno_nosurv,pat_GR,mut)

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({
for (rngNum in 1:30) {
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))
	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
											  setSeed=rngNum*5)
	# --------------------------------------------
	# feature selection - train only
	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
	dats_train <- lapply(dats,function(x) { 
						 x[,which(colnames(x) %in% pheno$ID)]})
	pat_GR_train <- pat_GR_all[which(pat_GR_all$ID %in% pheno$ID)]
	
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
	names(path_GRList) <- paste("MUT_",names(path_GRList),sep="")
	
	# group by pathway
	netList <- makePSN_NamedMatrix(dats_train$rna, rownames(dats_train$rna),
								   pathwayList,netDir,verbose=FALSE, 
								   numCores=numCores,writeProfiles=TRUE) 
	cat("Made RNA pathway nets\n")
	
	# each clinical var is its own net
	clinList <- list(age="age",gender="gender")
	netList2 <- makePSN_NamedMatrix(dats_train$clinical, 
									rownames(dats_train$clinical),
			clinList,netDir, simMetric="custom",customFunc=normDiff,
			sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
	cat("Made clinical nets\n")

	# add somatic mutations at pathway-level
	netList3 <- makePSN_RangeSets(pat_GR_train, path_GRList, netDir,
		numCores=numCores)
	cat("Made somatic mutation pathway nets\n")

	netList <- unlist(c(netList,netList2,netList3)) 
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
		tmp <- makePSN_NamedMatrix(dats$rna, rownames(dats$rna),
	        pathwayList[which(names(pathwayList)%in%pTally)],
			netDir,verbose=F,numCores=numCores)
	
		netList2 <- makePSN_NamedMatrix(dats$clinical, rownames(dats$clinical),
			clinList,netDir, simMetric="custom",customFunc=normDiff,
			sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)

		# add somatic mutations at pathway-level
		netList3 <- makePSN_RangeSets(pat_GR_all, path_GRList, netDir,
			numCores=numCores)
	
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

}}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})


