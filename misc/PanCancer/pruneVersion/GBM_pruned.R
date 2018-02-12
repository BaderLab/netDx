#' PanCancer binarized survival: GBM: Feature selection with one net per
#' datatype
#' 10-fold CV predictor design 
#' multi cutoff evaluation
#' also pruning RNA before running

rm(list=ls())
require(netDx)
require(netDx.examples)
source("../runLM.R")

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8

rootDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM"
inDir <- sprintf("%s/input",rootDir)
outRoot <- sprintf("%s/output",rootDir)

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/prune_%s",outRoot,dt)

# ----------------------------------------------------------------
# helper functions
# takes average of normdiff of each row in x
normDiff2 <- function(x) {
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

	sim <- matrix(0,nrow=ncol(x),ncol=ncol(x))
	for (k in 1:nrow(x)) {
		tmp <- normDiff(x[k,,drop=FALSE])
		sim <- sim + tmp
		rownames(sim) <- rownames(tmp)
		colnames(sim) <- colnames(tmp)
	}
	sim <- sim/nrow(x)
	sim
}

# -----------------------------------------------------------
# process input
inFiles <- list(
	clinical=sprintf("%s/GBM_clinical_core.txt",inDir),
	survival=sprintf("%s/GBM_binary_survival.txt",inDir)
	)
datFiles <- list(
	rna=sprintf("%s/GBM_mRNA_core.txt",inDir),
	mir=sprintf("%s/GBM_miRNA_core.txt",inDir),
	dnam=sprintf("%s/GBM_methylation_core.txt",inDir),
	cnv=sprintf("%s/GBM_CNV_core.txt",inDir)
)

pheno <- read.delim(inFiles$clinical,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- "ID"
# ------------------

surv <- read.delim(inFiles$survival,sep="\t",h=T,as.is=T)
colnames(surv)[1:2] <- c("ID","STATUS_INT")
survStr <- rep(NA,nrow(surv))
survStr[surv$STATUS_INT<1] <- "SURVIVENO"
survStr[surv$STATUS_INT>0] <- "SURVIVEYES"
surv$STATUS <- survStr
pheno <- merge(x=pheno,y=surv,by="ID")
pheno$X <- NULL
pheno_nosurv <- pheno[1:4]

cat("Collecting patient data:\n")
dats <- list() #input data in different slots
cat("\t* Clinical\n")
clinical <- pheno_nosurv
rownames(clinical) <- clinical[,1];
# =======================
# GBM-specific variables
clinical$performance_score[which(clinical$performance_score == "[Not Available]")] <- NA
clinical$performance_score <- strtoi(clinical$performance_score)
clinical$gender <- ifelse(pheno$gender=="FEMALE",1, 0)
# =======================
clinical$ID <- NULL
clinical <- t(clinical)
dats$clinical <- clinical; rm(clinical)

# create master input net
for (nm in names(datFiles)) {
	cat(sprintf("\t* %s\n",nm))
	tmp <- read.delim(datFiles[[nm]],sep="\t",h=T,as.is=T)
	if (colnames(tmp)[ncol(tmp)]=="X") tmp <- tmp[,-ncol(tmp)]
	rownames(tmp) <- tmp[,1]
	tmp <- t(tmp[,-1])
	class(tmp) <- "numeric"
	dats[[nm]] <- tmp
}

cat("\t Ordering column names\n")
# include only data for patients in classifier
dats <- lapply(dats, function(x) { x[,which(colnames(x)%in%pheno$ID)]})
dats <- lapply(dats, function(x) { 
	midx <- match(pheno$ID,colnames(x))
	x <- x[,midx]
	x
})

# confirm patient order the same for all input nets
pname <- colnames(dats[[1]])
for (k in 2:length(dats)) {
	if (all.equal(colnames(dats[[k]]),pname)!=TRUE) {
		cat(sprintf("Patient order doesn't match for %s\n",
			names(dats)[k]))
		browser()
	} 
}



# input nets for each category
netSets <- lapply(dats, function(x) rownames(x)) 

# compile data
alldat <- do.call("rbind",dats)
pheno_all <- pheno

combList <- list(    
    clinicalAcnv=c("clinical_cont","cnv.profile"),    
    clinical="clinical_cont",    
	mir="mir.profile",
	rna="rna.profile",
	cnv="cnv.profile",
	dnam="dnam.profile",
    clinicalArna=c("clinical_cont","rna.profile"),
    clinicalAmir=c("clinical_cont","mir.profile"),    
    clinicalAdnam=c("clinical_cont","dnam.profile"),    
    all="all"  
)

cat(sprintf("Clinical variables are: { %s }\n", 
	paste(rownames(dats$clinical),sep=",",collapse=",")))
rm(pheno,pheno_nosurv)

# ----------------------------------------------------------
# build classifier
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({

#### -----------------------------------------------------
### BEGIN PRUNING CODE
# apply pruning to proteomic data 
curwd <- getwd()
setwd("..")
source("LMprune.R")
source("runLM.R")
source("silh.R")
require(cluster)
setwd(curwd)
for (nm in setdiff(names(dats),"clinical")) {
print(nm)
	if (nrow(dats[[nm]])>10000) topVar <- 50 else topVar <- 100
	pdf(sprintf("%s/%s_prune.pdf",megaDir,nm))
	prune <- LMprune(dats[[nm]],pheno_all$STATUS,topVar=topVar)
	dev.off()
	if (!is.na(prune)) {
		if (prune$bestThresh < 1) {
		res <- prune$res
		res <- subset(res, adj.P.Val < prune$bestThresh)
		tmp <- dats[[nm]];orig_ct <- nrow(tmp)
		tmp <- tmp[which(rownames(tmp)%in% rownames(res)),]
		dats[[nm]] <- tmp
		netSets[[nm]] <- rownames(tmp)
		cat(sprintf("%s: Pruning with cutoff %1.2f\n", nm,prune$bestThresh))
		cat(sprintf("\t%i of %i left\n", nrow(tmp),orig_ct))
		}
	} else {
		cat(sprintf("%s: not pruning\n",nm))
	}
}
browser()
#### ----------------------------------------------------------


## Create the mega database with all patients and all nets.
## This will be used to predict test samples by subsetting just for feature
## selected nets in a given round
## Note that this is useful for all train/test splits because we can always
## change which samples are query and can always subset based on which nets
## are feature selected in a given round.
netDir <- sprintf("%s/networks",megaDir)
nonclin <- setdiff(names(netSets),"clinical")
netList <- makePSN_NamedMatrix(alldat,
	rownames(alldat),netSets[nonclin],netDir,
	verbose=FALSE,numCores=numCores,writeProfiles=TRUE)
netList2 <- makePSN_NamedMatrix(alldat, 
	rownames(alldat),netSets["clinical"],
	netDir,simMetric="custom",customFunc=normDiff2,
	verbose=FALSE,numCores=numCores,
	sparsify=TRUE,append=TRUE)
netList <- c(netList,netList2)
cat(sprintf("Total of %i nets\n", length(netList)))
browser()
	
# now create database
megadbDir	<- GM_createDB(netDir, pheno_all$ID, megaDir,numCores=numCores)

# first loop - over train/test splits
for (rngNum in 1:100) {
	rng_t0 <- Sys.time()
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))
	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
											  setSeed=rngNum*5)
	write.table(pheno_all,file=sprintf("%s/tt_split.txt",outDir),sep="\t",
		col=T,row=F,quote=F)
	# --------------------------------------------
	# feature selection - train only
	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
	alldat_train <- alldat[,which(colnames(alldat) %in% pheno$ID)]
	
	netDir <- sprintf("%s/networks",outDir)
	nonclin <- setdiff(names(netSets),"clinical")
	netList <- makePSN_NamedMatrix(alldat_train, 
		rownames(alldat_train),netSets[nonclin],
		netDir,verbose=FALSE,numCores=numCores,
		writeProfiles=TRUE)
	netList2 <- makePSN_NamedMatrix(alldat_train, 
		rownames(alldat_train),netSets["clinical"],
		netDir,simMetric="custom",customFunc=normDiff2,
		verbose=FALSE,numCores=numCores,
		sparsify=TRUE,append=TRUE)
	netList <- c(netList,netList2)
	cat(sprintf("Total of %i nets\n", length(netList)))
	
	# now create database
	dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)

	# second loop - over combinations of input data
 	for (cur in  names(combList)) {
		t0 <- Sys.time()
	    cat(sprintf("%s\n",cur))
	    pDir <- sprintf("%s/%s",outDir, cur)
	    dir.create(pDir)
	
		# run featsel once per subtype
		subtypes <- unique(pheno$STATUS)
		# run 10-fold cv per subtype
		for (g in subtypes) {
		    pDir2 <- sprintf("%s/%s",pDir,g)
		    if (file.exists(pDir2)) unlink(pDir2,recursive=TRUE)
			dir.create(pDir2)
		
			cat(sprintf("\n******\nSubtype %s\n",g))
			pheno_subtype <- pheno
			## label patients not in the current class as residual
			nong <- which(!pheno_subtype$STATUS %in% g)
			pheno_subtype$STATUS[nong] <- "nonpred"
			## sanity check
			print(table(pheno_subtype$STATUS,useNA="always"))
			resDir    <- sprintf("%s/GM_results",pDir2)
			## query for feature selection comprises of training 
			## samples from the class of interest
			trainPred <- pheno_subtype$ID[
				which(pheno_subtype$STATUS %in% g)]
			
			# Cross validation
			GM_runCV_featureSet(trainPred, resDir, dbDir$dbDir, 
				nrow(pheno_subtype),incNets=combList[[cur]],
				verbose=T, numCores=numCores,
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

		for (cutoff in 7:9) {
			predRes <- list()
			for (g in subtypes) {
				pDir2 <- sprintf("%s/%s",pDir,g)
				# get feature selected net names
				pTally <- read.delim(
				  sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir2,g),
				  sep="\t",h=T,as.is=T)
	
				# feature selected nets pass cutoff threshold
				pTally <- pTally[which(pTally[,2]>=cutoff),1]
				cat(sprintf("%s: %i pathways\n",g,length(pTally)))
				if (length(pTally)>=1) {
				curD <- sprintf("%s/cutoff%i",pDir2,cutoff)
				dir.create(curD)
				# query of all training samples for this class
				qSamps <- pheno_all$ID[which(pheno_all$STATUS %in% g & 
										 pheno_all$TT_STATUS%in%"TRAIN")]
			
				qFile <- sprintf("%s/%s_query",curD,g)
				GM_writeQueryFile(qSamps,incNets=pTally,
					nrow(pheno_all),qFile)
				resFile <- runGeneMANIA(megadbDir$dbDir,qFile,resDir=curD)
				predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),
					pheno_all,g)
				} else {
					predRes[[g]] <- NA
				}
			}
		
		oD <- sprintf("%s/cutoff%i",pDir,cutoff)
		dir.create(oD)
		outFile <- sprintf("%s/predictionResults.txt",oD)
		if (any(is.na(predRes))) {
			cat("One or more groups had zero feature selected nets\n")
			cat("# no feature-selected nets.\n",file=outFile) 
		}else {
			predClass <- GM_OneVAll_getClass(predRes)
			out <- merge(x=pheno_all,y=predClass,by="ID")
			write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)
			
			acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
			cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
				nrow(out), acc*100))
			
			require(ROCR)
			ROCR_pred <- prediction(out$SURVIVEYES_SCORE-out$SURVIVENO,
								out$STATUS=="SURVIVEYES")
			save(predRes,ROCR_pred,file=sprintf("%s/predRes.Rdata",oD))
		}
		}
		}
        
    #cleanup to save disk space
    system(sprintf("rm -r %s/dataset %s/tmp %s/networks", 
        outDir,outDir,outDir))
    system(sprintf("rm -r %s/dataset %s/networks", 
        outDir,outDir))

}
	pheno_all$TT_STATUS <- NA
	rng_t1 <- Sys.time()
	cat(sprintf("Time for one train/test split:"))
	print(rng_t1-rng_t0)

}, error=function(ex){
	print(ex)
}, finally={
	sink(NULL)
})
