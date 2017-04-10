#'  generate consensus profiles
rm(list=ls())
require(netDx)
require(netDx.examples)
require(ROCR)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9

inDir <- "/mnt/data2/BaderLab/PanCancer_LUSC/input"
outRoot <-"/mnt/data2/BaderLab/PanCancer_LUSC/output"
#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_LUSC/input"
#outRoot <-"/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_LUSC/output"

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/consNets_%s",outRoot,dt)
consNetDir <- "/mnt/data2/BaderLab/PanCancer_common"

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
	rppa=sprintf("%s/LUSC_RPPA_core.txt",inDir),
	mut=sprintf("%s/from_firehose/LUSC_core_somatic_mutations.txt",
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
dats$clinical <- clin; rm(clin)

#### change this for current tumour type
clinList <- list(age="age",stage="stage")

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

# RNA
cat("\t* RNA\n")
rna <- read.delim(inFiles$rna,sep="\t",h=T,as.is=T)
rna <- t(rna)
colnames(rna) <- rna[1,]; rna <- rna[-1,]; 
rna <- rna[-nrow(rna),]
class(rna) <- "numeric"
rownames(rna) <- sub("mRNA_","",rownames(rna))
dpos <- regexpr("\\.", rownames(rna))
rownames(rna) <- substr(rownames(rna),1,dpos-1)
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

rm(pheno,pat_GR,mut)

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({
#for (rngNum in 4) { #1:100) {
###	cat(sprintf("-------------------------------\n"))
###	cat(sprintf("RNG seed = %i\n", rngNum))
###	cat(sprintf("-------------------------------\n"))
###	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
###	dir.create(outDir)

outDir <- megaDir
	pheno_all$TT_STATUS <- "DUMMY"
###	# --------------------------------------------
###	# feature selection - train only
###	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
###	dats_train <- lapply(dats,function(x) { 
###						 x[,which(colnames(x) %in% pheno$ID)]})
###	pat_GR_train <- pat_GR_all[which(pat_GR_all$ID %in% pheno$ID)]
###	
	# create nets 
	netDir <- sprintf("%s/networks",outDir) 
	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
	   path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)
	PROT_pathwayList <- pathwayList
	names(PROT_pathwayList) <- paste("PROT", names(pathwayList),sep="_")

	data(genes)
	gene_GR     <- GRanges(genes$chrom,IRanges(genes$txStart,genes$txEnd),
	  	name=genes$name2)
	cat("* Limiting to pathway genes\n")
	path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
	names(path_GRList) <- paste("MUT_",names(path_GRList),sep="")
	
###	# RNA - group by pathway
###	netList <- makePSN_NamedMatrix(dats_train$rna, rownames(dats_train$rna),
###								   pathwayList,netDir,verbose=FALSE, 
###								   numCores=numCores,writeProfiles=TRUE) 
###	cat("Made RNA pathway nets\n")
###
###	# PROTEIN - group by pathway
###	netList2 <- makePSN_NamedMatrix(dats_train$rppa, rownames(dats_train$rppa),
###   	     PROT_pathwayList,netDir,verbose=FALSE, 
###   	     numCores=numCores,writeProfiles=TRUE,append=TRUE) 
###	cat("Made protein pathway nets\n")
###	
###	# each clinical var is its own net
###	netList3 <- makePSN_NamedMatrix(dats_train$clinical, 
###									rownames(dats_train$clinical),
###			clinList,netDir, simMetric="custom",customFunc=normDiff,
###			sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
###	cat("Made clinical nets\n")
###
###	# add somatic mutations at pathway-level
###	netList4 <- makePSN_RangeSets(pat_GR_train, path_GRList, netDir,
###		numCores=numCores)
###	cat("Made somatic mutation pathway nets\n")
###
###	netList <- unlist(c(netList,netList2,netList3,netList4)) 
###	cat(sprintf("Total of %i nets\n", length(netList)))
###	
###	# now create database
###	dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)
###	
	# run featsel once per subtype
	subtypes <- unique(pheno_all$STATUS)
	# run 10-fold cv per subtype
###	for (g in subtypes) {
###	    pDir <- sprintf("%s/%s",outDir,g)
###	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
###		dir.create(pDir)
###	
###		cat(sprintf("\n******\nSubtype %s\n",g))
###		pheno_subtype <- pheno
###		## label patients not in the current class as a residual
###		pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
###		## sanity check
###		print(table(pheno_subtype$STATUS,useNA="always"))
###		resDir    <- sprintf("%s/GM_results",pDir)
###		## query for feature selection comprises of training 
###		## samples from the class of interest
###		trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% g)]
###		
###		# Cross validation
###		GM_runCV_featureSet(trainPred, resDir, dbDir$dbDir, 
###			nrow(pheno_subtype),verbose=T, numCores=numCores,
###			GMmemory=GMmemory)
###	
###		# patient similarity ranks
###		prank <- dir(path=resDir,pattern="PRANK$")
###		# network ranks
###		nrank <- dir(path=resDir,pattern="NRANK$")
###		cat(sprintf("Got %i prank files\n",length(prank)))
###			
###	    # Compute network score
###		pTally		<- GM_networkTally(paste(resDir,nrank,sep="/"))
###		head(pTally)
###		# write to file
###		tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
###		write.table(pTally,file=tallyFile,sep="\t",col=T,row=F,quote=F)
###	}
	
	## ----class-prediction, eval=TRUE-------------------------
	# now create GM databases for each class
	# should contain train + test patients
	# and be limited to nets that pass feature selection
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		if (!file.exists(pDir)) dir.create(pDir)

		# get feature selected net names
		pTally <- read.delim(
			sprintf("%s/LUSC_%s_consNets.txt",consNetDir,g),
			sep="\t",h=T,as.is=T)[,1]
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
	   		     pathwayList[idx],writeProfiles=TRUE,
				netDir,verbose=F,numCores=numCores)
		}
	
		# clinical
		idx <- which(names(clinList) %in% pTally)
		if (any(idx)) {
			cat(sprintf("clinical: included %i nets\n", length(idx)))
		netList2 <- makePSN_NamedMatrix(dats$clinical, rownames(dats$clinical),
			clinList[idx],
			netDir, simMetric="custom",customFunc=normDiff,
			sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
		}

		# add somatic mutations at pathway-level
		idx <- which(names(path_GRList) %in% pTally) 
		if (any(idx)) {
			cat(sprintf("mutations: included %i nets\n", length(idx)))
			netList3 <- makePSN_RangeSets(pat_GR_all, 
				path_GRList[idx], 
				netDir,numCores=numCores)
		}

		# proteomics group by pathway
		idx <- which(names(PROT_pathwayList) %in% pTally)
		if (any(idx)) {
			cat(sprintf("proteomics: included %i nets\n", length(idx)))
			netList4 <- makePSN_NamedMatrix(dats$rppa, 
				rownames(dats$rppa),
   		     	PROT_pathwayList[idx],netDir,verbose=FALSE, 
   		     	numCores=numCores,writeProfiles=TRUE,append=TRUE) 
		}

###		# create db
###		dbDir <- GM_createDB(netDir,pheno$ID,pDir,numCores=numCores)
###		# query of all training samples for this class
###		qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
###								 pheno$TT_STATUS%in%"TRAIN")]
###	
###		qFile <- sprintf("%s/%s_query",pDir,g)
###		GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
###		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
###		predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
	}
	
###	predClass <- GM_OneVAll_getClass(predRes)
###	out <- merge(x=pheno_all,y=predClass,by="ID")
###	outFile <- sprintf("%s/predictionResults.txt",outDir)
###	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
###	
###	acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
###	cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
###		nrow(out), acc*100))
###	
###	ROCR_pred <- prediction(out$SURVIVEYES_SCORE-out$SURVIVENO,
###						out$STATUS=="SURVIVEYES")
###	save(predRes,ROCR_pred,file=sprintf("%s/predRes.Rdata",outDir))

###	# cleanup
###	system(sprintf("rm -r %s/dataset %s/tmp %s/networks", 
###		outDir,outDir,outDir))	
###	system(sprintf("rm -r %s/SURVIVENO/dataset %s/SURVIVENO/networks",
###		outDir,outDir))
###	system(sprintf("rm -r %s/SURVIVEYES/dataset %s/SURVIVEYES/networks",
###		outDir,outDir))
###	system(sprintf("rm -r %s/SURVIVEYES/tmp %s/SURVIVENO/tmp",
###		outDir,outDir))

#}
}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})


