#' PanCancer predictor: univariate filtering by lasso + gene-level nets
#' similarity by Euclidean distance + local scaling

# ----------------------------------------------------------------
# helper functions
sim.pearscale <- function (dat, K = 20, alpha = 0.5) {
ztrans <- function(m) {
	m <- as.matrix(m)
	m2 <- apply(m,1,function(x) { (x-mean(x,na.rm=T))/sd(x,na.rm=T)})
	m2
}
normalize <- function(X) {
		print(dim(X))
        row.sum.mdiag <- rowSums(X,na.rm=T) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
		X <- (X+t(X))/2
        return(X)
}

	if (nrow(dat)<6) {	
		z1 <- ztrans(dat)
		euc <- as.matrix(dist(z1,method="euclidean"))^(1/2)
	} else {
		euc <- as.matrix(1-cor(dat,method="pearson",use="pairwise.complete.obs"))
	}
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort,na.last=TRUE)))
	print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)],na.rm=T))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) +
        .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)

    W <- (densities + t(densities))/2
	W <- normalize(W)
	# remove patients with no datapoints (full column/row of NAs)
	idx <- which(rowSums(is.na(euc))==ncol(W)-1)
	if (any(idx)) { 
		W <- W[-idx,]
		idx <- which(colSums(is.na(euc))==ncol(W)-1)
		W <- W[,-idx]
	}
    return(W)
}

# ----------------------------------------------------------------
runPredictor <- function(mega_combList,rngVals,netSets,dats,pheno_all,megaDir,
	cutoffSet) {
require(netDx)
require(netDx.examples)
require(glmnet)

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

	# lasso
	for (nm in names(dats_train)) { 
		print(nm)
		if (nrow(dats_train[[nm]])<2)  # clinical only has one var, take it.
			vars <- rownames(dats_train[[nm]])
		else { 
			fit <- cv.glmnet(x=t(dats_train[[nm]]),
				y=factor(pheno$STATUS), family="binomial", alpha=1) # lasso
			wt <- abs(coef(fit,s="lambda.min")[,1])
			vars <- setdiff(names(wt)[which(wt>.Machine$double.eps)],"(Intercept)")
		}
		cat(sprintf("rngNum %i: %s: %s pruned\n",rngNum,nm,length(vars)))
		if (length(vars)>0) {
			tmp <- dats_train[[nm]]
			tmp <- tmp[which(rownames(tmp) %in% vars),,drop=FALSE]
			dats_train[[nm]] <- tmp
			for (k in rownames(tmp)) { netSets_iter[[k]] <- k }
			combList[[nm]] <- sprintf("%s_cont", rownames(tmp))
		} else {
			# leave dats_train as is, make a single net
			netSets_iter[[nm]] <- rownames(dats_train[[nm]])
			combList[[nm]] <- sprintf("%s_cont",nm)
		} 
	}
	
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

	# END lasso UF
	# ----------------------
	alldat_train <- do.call("rbind",dats_train)
	netDir <- sprintf("%s/networks",outDir)

	cat(sprintf("Making test nets for rng%i\n", rngNum))
	netList <- makePSN_NamedMatrix(alldat_train,
        rownames(alldat_train),netSets_iter,netDir,
        simMetric="custom",customFunc=sim.pearscale,
        writeProfiles=FALSE,
        sparsify=TRUE,useSparsify2=TRUE,cutoff=.Machine$double.eps,
		sparsify_edgeMax=maxEdge,
        verbose=FALSE,numCores=numCores)
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
	netList <- makePSN_NamedMatrix(alldat_tmp,
        rownames(alldat_tmp),netSets_iter,netDir,
        simMetric="custom",customFunc=sim.pearscale,
        writeProfiles=FALSE,
        sparsify=TRUE,useSparsify2=TRUE,cutoff=.Machine$double.eps,
		sparsify_edgeMax=maxEdge,
        verbose=TRUE,numCores=numCores)
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
