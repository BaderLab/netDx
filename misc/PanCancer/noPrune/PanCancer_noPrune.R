#' PanCancer binarized survival: KIRC: Feature selection with one net per
# datatype
#' 10-fold CV predictor design 

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

# ----------------------------------------------------------------
runPredictor <- function(mega_combList,rngVals,netSets,dats,pheno_all,megaDir,
	cutoffSet,maxEdge,spCutoff) {

require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- cutoffSet
cat(sprintf("FS cutoff = %i\n", cutoff))

if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({

# first loop - over train/test splits
for (rngNum in rngVals) {
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

	dats_train <- lapply(dats, function(x) x[,which(colnames(x) %in% pheno$ID),
			drop=F])
	alldat_train <- do.call("rbind",dats_train)
	netSets_iter <- netSets
	
	netDir <- sprintf("%s/networks",outDir)
	cat(sprintf("Making test nets for rng%i\n", rngNum))
	netLen <- unlist(lapply(netSets_iter,length))
	pearnet <- which(netLen>=6)
	othernet <- setdiff(1:length(netSets_iter),pearnet)
	netList <- c(); netList2 <- c()
	if (any(othernet)) {
	netList <- makePSN_NamedMatrix(alldat_train,
        rownames(alldat_train),netSets_iter[othernet],netDir,
        simMetric="custom",customFunc=normDiff2,
        writeProfiles=FALSE,
        sparsify=TRUE,useSparsify2=TRUE,cutoff=spCutoff,
		sparsify_edgeMax=maxEdge,
        verbose=FALSE,numCores=numCores)
	}
	if (any(pearnet)) {
		netList2 <- makePSN_NamedMatrix(alldat_train,
	        rownames(alldat_train),netSets_iter[pearnet],netDir,
	        writeProfiles=TRUE,
	        verbose=FALSE,numCores=numCores,append=TRUE)
	}
	netList <- c(netList,netList2)
	cat(sprintf("Total of %i nets\n", length(netList)))

	# now create database
	dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)

	# second loop - over combinations of input data
 	for (cur in names(combList)) {
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

	netDir <- sprintf("%s/test_networks",outDir)
	netLen <- unlist(lapply(netSets_iter,length))
	pearnet <- which(netLen>=6)
	othernet <- setdiff(1:length(netSets_iter),pearnet)
	netList <- c(); netList2 <- c()
	alldat <- do.call("rbind",dats)
	if (any(othernet)) {
	netList <- makePSN_NamedMatrix(alldat,
        rownames(alldat),netSets_iter[othernet],netDir,
        simMetric="custom",customFunc=normDiff2,
        writeProfiles=FALSE,
        sparsify=TRUE,useSparsify2=TRUE,cutoff=spCutoff,
		sparsify_edgeMax=maxEdge,
        verbose=FALSE,numCores=numCores)
	}
	if (any(pearnet)) {
		netList2 <- makePSN_NamedMatrix(alldat,
	        rownames(alldat),netSets_iter[pearnet],netDir,
	        writeProfiles=TRUE,
	        verbose=FALSE,numCores=numCores,append=TRUE)
	}
	netList <- c(netList,netList2)
	cat(sprintf("Total of %i nets\n", length(netList)))

	# now create database
	testdbDir	<- GM_createDB(netDir, pheno_all$ID, 
		outDir,numCores=numCores)

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

			# query of all training samples for this class
			qSamps <- pheno_all$ID[which(pheno_all$STATUS %in% g & 
				 pheno_all$TT_STATUS%in%"TRAIN")]
		
			qFile <- sprintf("%s/%s_query",pDir2,g)
			GM_writeQueryFile(qSamps,incNets=pTally,
				nrow(pheno_all),qFile)
			resFile <- runGeneMANIA(testdbDir$dbDir,qFile,resDir=pDir2)
			predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),
				pheno_all,g)
			} else {
				predRes[[g]] <- NA
			}
		}
		
		if (any(is.na(predRes))) {
			cat("One or more groups had zero feature selected nets\n")
			cat("# no feature-selected nets.\n",file=outFile) 
		}else {
			predClass <- GM_OneVAll_getClass(predRes)
			out <- merge(x=pheno_all,y=predClass,by="ID")
			outFile <- sprintf("%s/predictionResults.txt",pDir)
			write.table(out,file=outFile,sep="\t",col=T,row=F,
				quote=F)
			
			acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
			cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
				nrow(out), acc*100))
			
			require(ROCR)
			ROCR_pred <- prediction(out$SURVIVEYES_SCORE-out$SURVIVENO,
					out$STATUS=="SURVIVEYES")
			save(predRes,ROCR_pred,file=sprintf("%s/predRes.Rdata",pDir))
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
}
