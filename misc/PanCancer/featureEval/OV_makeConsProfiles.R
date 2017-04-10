#' feature selection for LUSC from PanCancer survival dataset
#' 10-fold CV predictor design with clinical and mRNA data
rm(list=ls())
require(netDx)
require(netDx.examples)
require(ROCR)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9
maxRng <- 25 

inDir <- "/mnt/data2/BaderLab/PanCancer_OV/input"
outRoot <-"/mnt/data2/BaderLab/PanCancer_OV/output"

#### makeConsProfiles code block >>>>
consNetDir <- "/mnt/data2/BaderLab/PanCancer_common"
maxRng <- 25 

analysisMode <- "consNets" # none |consNets | bestConsNets | randomNets
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

### >>> CHANGE THIS for the current tumour
clinList <- list(age="age") # clinical nets
protList <- list(NOTCH="NOTCH3") # proteomics nets
# probe-> gene mappings from BioC package IlluminaHumanMethylation27k.db,
# as.list(IlluminaHumanMethylation27kALIAS2PROBE) function
methList <- list(BRCA1=c("cg06973652","cg11964474"),
				 BRCA2=c("cg12836863","cg27253386"))

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
	clinical=sprintf("%s/OV_clinical_core.txt",inDir),
	survival=sprintf("%s/OV_binary_survival.txt",inDir),
	rna=sprintf("%s/OV_mRNA_core.txt",inDir),
	rppa=sprintf("%s/OV_RPPA_core.txt",inDir),
	mut=sprintf("%s/from_firehose/OV_core_somatic_mutations.txt",
		inDir),
	meth=sprintf("%s/OV_methylation_short.txt",inDir) # pre-prepared short file
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
rownames(clin) <- clin[,1]; 
clin <- t(clin[,2,drop=FALSE])
dats$clinical <- clin; rm(clin)

# methylation
cat("\t* DNA methylation\n")
meth <- read.delim(inFiles$meth,sep="\t",h=T,as.is=T)
dats$meth <- meth; rm(meth)
colnames(dats$meth) <- gsub("\\.","-",colnames(dats$meth))

### done once, don't repeat.
###rownames(meth) <- meth[,1]; meth <- meth[,-1]
###meth <- t(meth)
###rownames(meth) <- sub("methylation_","",rownames(meth))
###meth <- meth[which(rownames(meth) %in% unlist(methList)),]
###write.table(meth, file=sprintf("%s/OV_methylation_short.txt", inDir),sep="\t",
###	col=TRUE,row=TRUE)

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
dats$rna <- rna; rm(rna)

# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID),
	drop=FALSE]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx,drop=FALSE]
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
	outDir <- megaDir 
	pheno_all$TT_STATUS <- "DUMMY"
	
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
		if (!file.exists(pDir)) dir.create(pDir)
#### makeConsProfiles code block >>>>
		pTally <- switch(analysisMode,
	       consNets={
			read.delim(
				sprintf("%s/OV_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
		}, none={
			read.delim(
				sprintf("%s/OV_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
		}, bestConsNets={
			tmp <- read.delim(
				sprintf("%s/OV_%s_consNets.txt", consNetDir,g),
				sep="\t",h=T,as.is=T)[,1]
			fName <- dir(consNetDir,pattern=sprintf("OV_%s_correlations",g))
			bestSet <- read.delim(
				sprintf("%s/%s",consNetDir,fName),sep="\t",h=T,as.is=T)
			idx_tmp <- unlist(sapply(1:nrow(bestSet), 
				function(k) max(bestSet[k,4:6])))
			oldlen <- length(tmp)
			tmp <- tmp[which(tmp %in% rownames(bestSet)[which(idx_tmp>2)])]
			cat(sprintf("Filtered %i to %i best\n", oldlen,length(tmp)))
			tmp
		}, randomNets= {
			cat("need to count num nets in cons, sigh")
			numCons <- sprintf("wc -l %s/OV_%s_consNets.txt",consNetDir,g) 
			numCons <- as.integer(strsplit(system(numCons,intern=T)," ")[[1]][1])
			pTally <- read.delim(
				sprintf("%s/OV_%s_AllNets.txt", consNetDir,g),
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
		# RNA
        idx <- which(names(pathwayList) %in% pTally)
        if (any(idx)) {
            cat(sprintf("RNA: included %i nets\n", length(idx)))
            tmp <- makePSN_NamedMatrix(dats$rna, rownames(dats$rna),
                 pathwayList[idx],
                netDir,verbose=F,numCores=numCores, writeProfiles=TRUE)
        }
	
            
        idx <- which(names(clinList) %in% pTally)
		if (any(idx)) {
            cat(sprintf("clinical: included %i nets\n", length(idx)))
            netList2 <- makePSN_NamedMatrix(dats$clinical, rownames(dats$clinical),
                clinList[idx],
                netDir, simMetric="custom",customFunc=normDiff,
                sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
        }

		# somatic mutations - add somatic mutations at pathway-level
		idx <- which(names(path_GRList) %in% pTally) 
		if (any(idx)) {
			cat(sprintf("mutations: included %i nets\n", length(idx)))
			netList3 <- makePSN_RangeSets(pat_GR_all, 
				path_GRList[idx], 
				netDir,numCores=numCores)
        }
    
		# rppa by genes of interest  group by pathway
        idx <- which(names(protList) %in% pTally) 
        if (any(idx)){
            cat(sprintf("protein: included %i nets\n", length(idx)))
            netList4 <- makePSN_NamedMatrix(dats$rppa, 
                rownames(dats$rppa),
                 protList[idx],netDir,simMetric="custom", customFunc=normDiff,
                sparsify=TRUE,verbose=TRUE,numCores=numCores,append=TRUE)
        }
	
		# methylation BRCA1/2 only
        idx <- which(names(methList) %in% pTally) 
        if (any(idx)){
            cat(sprintf("methylation: included %i nets\n", length(idx)))
            netList5 <- makePSN_NamedMatrix(dats$meth, 
                rownames(dats$meth),methList[idx],netDir,verbose=FALSE, 
                numCores=numCores,writeProfiles=TRUE,append=TRUE) 
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
