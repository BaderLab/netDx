#' PanCancer predictor: univariate filtering by lasso + gene-level nets
#' similarity by Euclidean distance + local scaling

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
	cutoffSet) {
require(netDx)
require(netDx.examples)
require(pamr)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
maxEdge <- 6000  ### max edge after sparsification

if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
tryCatch({

alldat <- do.call("rbind",dats)

for (rngNum in rngVals) {
	combList <- mega_combList # clean slate
	rng_t0 <- Sys.time()
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))

	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)
	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
				  setSeed=rngNum*5)
	write.table(pheno_all,file=sprintf("%s/tt_split.txt",outDir),
		sep="\t",col=T,row=F,quote=F)

	# feature selection - train only
	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")
	dats_train <- lapply(dats, function(x) x[,which(colnames(x) %in% pheno$ID),
		drop=FALSE])
	netSets_iter <- list()

	# impute
	dats_train <- lapply(dats_train, function(x) {
		missidx <- which(rowSums(is.na(x))>0) 
		for (i in missidx) {
			na_idx <- which(is.na(x[i,]))
			x[i,na_idx] <- median(x[i,],na.rm=TRUE) 
		}
		x
	})

	# pamr
	for (nm in setdiff(names(dats_train),"clinical")) { 
		print(nm)
		if (nrow(dats_train[[nm]])<2)  # clinical only has one var, take it.
			vars <- rownames(dats_train[[nm]])
		else { 
			set.seed(123); # reproducible
			tmp <- na.omit(dats_train[[nm]])
			data <- list(x=tmp,y=factor(pheno$STATUS),
				genenames=rownames(tmp),
				geneid=rownames(tmp))
			data.fit <- pamr.train(data)
			data.cv <- pamr.cv(data.fit, data)
			idx <- which.min(data.cv$error)
			thresh <- data.cv$threshold[idx]
			vars <- pamr.listgenes(data.fit,data,thresh,data.cv)[,1]
		}
		cat(sprintf("rngNum %i: %s: %s pruned\n",rngNum,nm,length(vars)))
		tmp <- dats_train[[nm]];orig_ct <- nrow(tmp)
		tmp <- tmp[which(rownames(tmp)%in% vars),]
		dats_train[[nm]] <- tmp
		cat(sprintf("\t%i:pamr:%s:%i of %i left\n",rngNum,nm,nrow(tmp),orig_ct))
		netSets_iter[[nm]] <- rownames(dats_train[[nm]])
	}
	netSets_iter[["clinical"]] <- rownames(dats_train[["clinical"]])
	
	if ("clinicalArna" %in% names(combList)) 
		combList[["clinicalArna"]] <- c(combList[["clinical"]],combList[["rna"]])
	if ("clinicalAmir" %in% names(combList)) 
		combList[["clinicalAmir"]] <- c(combList[["clinical"]],combList[["mir"]])
	if ("clinicalAcnv" %in% names(combList)) 
		combList[["clinicalAcnv"]] <- c(combList[["clinical"]],combList[["cnv"]])
	if ("clinicalAdnam" %in% names(combList)) 
		combList[["clinicalAdnam"]] <- c(combList[["clinical"]],combList[["dnam"]])
	if ("clinicalAprot" %in% names(combList)) 
		combList[["clinicalAprot"]] <- c(combList[["clinical"]],combList[["prot"]])

	# END pamr UF
	# ----------------------
	alldat_train <- do.call("rbind",dats_train)
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
        sparsify=TRUE,useSparsify2=TRUE,cutoff=.Machine$double.eps,
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
	
	dbDir	<- GM_createDB(netDir, pheno$ID, outDir,numCores=numCores)

	# second loop - over combinations of input data
 	for (cur in  names(combList)) {
		t0 <- Sys.time()
	    cat(sprintf("Input datatype\n%s\n",cur))
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
			print(table(pheno_subtype$STATUS,useNA="always")) # sanitycheck
			resDir    <- sprintf("%s/GM_results",pDir2)
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

	# make test db
	# impute
	train_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TRAIN")]
	test_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TEST")]
	dats_tmp <- lapply(dats, function(x) {
		missidx <- which(rowSums(is.na(x))>0) 
		train_idx <- which(colnames(x) %in% train_samp)
		test_idx <- which(colnames(x) %in% test_samp)
		for (i in missidx) {
			# impute train and test separately
			na_idx <- intersect(which(is.na(x[i,])),train_idx)
			na_idx1 <- na_idx
			x[i,na_idx] <- median(x[i,train_idx],na.rm=TRUE) 

			na_idx <- intersect(which(is.na(x[i,])),test_idx)
			na_idx2 <- na_idx
			x[i,na_idx] <- median(x[i,test_idx],na.rm=TRUE) 
		}
		x
	})
	alldat_tmp <- do.call("rbind",dats_tmp)
	netDir <- sprintf("%s/test_networks",outDir)

	netLen <- unlist(lapply(netSets_iter,length))
	pearnet <- which(netLen>=6)
	othernet <- setdiff(1:length(netSets_iter),pearnet)
	netList <- c(); netList2 <- c()
	if (any(othernet)) {
	netList <- makePSN_NamedMatrix(alldat_tmp,
        rownames(alldat_tmp),netSets_iter[othernet],netDir,
        simMetric="custom",customFunc=normDiff2,
        writeProfiles=FALSE,
        sparsify=TRUE,useSparsify2=TRUE,cutoff=.Machine$double.eps,
		sparsify_edgeMax=maxEdge,
        verbose=FALSE,numCores=numCores)
	}
	if (any(pearnet)) {
		netList2 <- makePSN_NamedMatrix(alldat_tmp,
	        rownames(alldat_tmp),netSets_iter[pearnet],netDir,
	        writeProfiles=TRUE,
	        verbose=FALSE,numCores=numCores,append=TRUE)
	}
	netList <- c(netList,netList2)
	cat(sprintf("Total of %i nets\n", length(netList)))

	# now create database
	testdbDir	<- GM_createDB(netDir, pheno_all$ID, outDir,numCores=numCores)

		# classify patients
		for (cutoff in cutoffSet) {
			predRes <- list()
			for (g in subtypes) {
				pDir2 <- sprintf("%s/%s",pDir,g)

				pTally <- read.delim(
				  sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir2,g),
				  sep="\t",h=T,as.is=T)
				pTally <- pTally[which(pTally[,2]>=cutoff),1]
				cat(sprintf("%s: %i pathways\n",g,length(pTally)))
				if (length(pTally)>=1) {
					curD <- sprintf("%s/cutoff%i",pDir2,cutoff)
					dir.create(curD)
					qSamps <- pheno_all$ID[which(pheno_all$STATUS %in% g & 
											 pheno_all$TT_STATUS%in%"TRAIN")]
				
					qFile <- sprintf("%s/%s_query",curD,g)
					GM_writeQueryFile(qSamps,incNets=pTally,
						nrow(pheno_all),qFile)
					resFile <- runGeneMANIA(testdbDir$dbDir,qFile,resDir=curD)
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
		}
		}

	} # input data combinations
        
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
