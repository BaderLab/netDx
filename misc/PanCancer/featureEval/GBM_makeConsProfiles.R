#' write consensus nets for feature evaluation

rm(list=ls())
require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8

inDir <- "/home/netdx/BaderLab/PanCancer_GBM/input"
outRoot <-"/home/netdx/BaderLab/PanCancer_GBM/output"
#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/input"
#outRoot <-"/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/output"

#### makeConsProfiles code block >>>>
consNetDir <- "/home/netdx/BaderLab/PanCancer_common"
maxRng <- 25 

for (analysisMode in c("consNets","bestConsNets","randomNets")) {
#analysisMode <- "consNets" # none |consNets | bestConsNets | randomNets
# none = just generate nets of consensus profiles, don't run predictions
# consNets = predictions with GMdb of all nets
# bestConsNets = GM db with nets that have corr < 0.01
# randomNets = GMdb with randomly sampled nets of size consNets 

if (!analysisMode %in% c("none","consNets","bestConsNets","randomNets")){
	cat("invalid method")
}

if (!file.exists(outRoot)) dir.create(outRoot)
dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/%s_%s",outRoot,analysisMode,dt)
#### <<<< makeConsProfiles code block 


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

clinList <- list(age="age",gender="gender",Karnofsky="performance_score")

rm(pheno,pheno_nosurv,pat_GR,mut)

# ----------------------------------------------------------
# build classifier
numCores <- 8L
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({
	outDir <- megaDir
	dir.create(outDir)

	pheno_all$TT_STATUS <- "DUMMY"
###	# --------------------------------------------
	
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
	
	
	# run featsel once per subtype
	subtypes <- unique(pheno_all$STATUS)
	
	## ----class-prediction, eval=TRUE-------------------------
	# now create GM databases for each class
	# should contain train + test patients
	# and be limited to nets that pass feature selection
	pheno <- pheno_all
	predRes <- list()
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		if (!file.exists(pDir)) dir.create(pDir,recursive=TRUE)
#### makeConsProfiles code block >>>>
		pTally <- switch(analysisMode,
	       consNets={
			read.delim(
				sprintf("%s/GBM_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
		}, none={
			read.delim(
				sprintf("%s/GBM_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
		}, bestConsNets={
			tmp <- read.delim(
				sprintf("%s/GBM_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
			fName <- dir(consNetDir,pattern=sprintf("GBM_%s_correlations",g))
			bestSet <- read.delim(
				sprintf("%s/%s",consNetDir,fName),sep="\t",h=T,as.is=T)
			idx_tmp <- unlist(sapply(1:nrow(bestSet), 
				function(k) max(bestSet[k,4:6])))
			oldlen <- length(tmp)
			tmp <- tmp[which(tmp %in% rownames(bestSet)[which(idx_tmp>2)])]
			cat(sprintf("Filtered %i to %i best\n", oldlen,length(tmp)))
			tmp
		}, randomNets= {
			numCons <- sprintf("wc -l %s/GBM_%s_consNets.txt",consNetDir,g) 
			numCons <- as.integer(strsplit(system(numCons,intern=T)," ")[[1]][1])
			pTally <- read.delim(
				sprintf("%s/GBM_%s_AllNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
			cat(sprintf("**** Sampling %i nets randomly ****\n",numCons))
			set.seed(249);
			pTally <- sample(pTally, numCons, replace=FALSE)
		},{
			stop(sprintf("Invalid analysisMode = %s\n", analysisMode))
		})
#### <<< makeConsProfiles code block 
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		cat(sprintf("%s: %i pathways\n",g,length(pTally)))
		netDir <- sprintf("%s/networks",pDir)
	
		# prepare nets for new db
		# RNA - 
		idx <- which(names(pathwayList) %in% pTally)
		if (any(idx)) {
		netList <- makePSN_NamedMatrix(dats$rna, rownames(dats$rna),
	        pathwayList[idx],writeProfiles=TRUE,
			netDir,verbose=F,numCores=numCores)
			cat(sprintf("TEST: made %i RNA nets\n", length(netList)))
		}
	
		# clinical	
		idx <- which(names(clinList) %in% pTally)
		if (any(idx)) {	
			netList2 <- makePSN_NamedMatrix(dats$clinical, 
				rownames(dats$clinical),
				clinList[idx],netDir, simMetric="custom",customFunc=normDiff,
				sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
			cat(sprintf("TEST: made %i clinical nets\n", length(netList2)))
		}

		# add somatic mutations at pathway-level
		idx <- which(names(path_GRList) %in% pTally)
		if (any(idx)) {
			netList3 <- makePSN_RangeSets(pat_GR_all, 
				path_GRList[idx], 
				netDir,numCores=numCores)
			cat(sprintf("TEST: made %i mutation nets\n", length(netList3)))
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

	if (!analysisMode %in% "none") {
	resDir <- sprintf("%s/predictions",outDir)
	dir.create(resDir)
	for (k in 1:length(predRes)) {
		predClass <- GM_OneVAll_getClass(predRes[[k]])
 		out <- merge(x=pheno_all,y=predClass,by="ID")
		outFile <- sprintf("%s/predictionResults_%i.txt",resDir,k)
	 	write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
	}
	}
### <<<< makeConsProfiles code block 
}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})


}
