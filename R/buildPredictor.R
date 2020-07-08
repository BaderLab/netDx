#' Run nested cross-validation on data 
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' 
#' @details wrapper function to run netDx with nested cross-validation, 
#' with an inner loop of X-fold cross-validation and an outer loop of different
#' random splits of data into train and blind test. The user needs to supply
#' a custom function to create PSN, see createPSN_MultiData(). This wrapper
#' provides flexibility for designs with one or several heterogeneous data
#' types, and one or more ways of defining patient similarity. 
#' For example, designs it handles includes
#' 1) Single datatype, single similarity metric: Expression data -> pathways
#' 2) Single datatype, multiple metrics: Expression data -> pathways
#'	(Pearson corr) and single gene networks (normalized difference)
#' 3) Multiple datatypes, multiple metrics: Expression -> Pathways; 
#'	Clinical -> single or grouped nets
#' @param dataList (MultiAssayExperiment) sample metadata. Clinical data is 
#' in colData() and other input datatypes are in assays() slot.
#' names(groupList) should match names(assays(dataList)). The only exception
#' is clinical data. If a groupList entry is called "clinical", the algorithm
#' will search for corresponding variable names in colData(dataList) (i.e.
#' columns of sample metadata table).
#' @param groupList (list of lists) keys are datatypes, and values are 
#' lists indicating how units for those datatypes are to be grouped. 
#' Keys must match names(assays(dataList)). The only exception is for clinical
#' values. Variables for "clinical" will be extracted from columns of the 
#' sample metadata table (i.e. from colData(dataList)).   
#' e.g. groupList[["rna"]] could be a list of pathway definitions. 
#' So keys(groupList[["rna"]]) would have pathway names, generating one PSN
#' per pathways, and values(groupList[["rna"]]) would be genes that would be
#' grouped for the corresponding pathwayList.
#' @param makeNetFunc (function) user-defined function for creating the set
#' of input PSN provided to netDx. See createPSN_MultiData()::customFunc.
#' @param outDir (char) directory where results will be stored. If this 
#' directory exists, its contents will be overwritten
#' @param trainProp (numeric 0 to 1) Percent samples to use for training
#' @param featScoreMax (integer) number of CV folds in inner loop
#' @param numSplits (integer) number of train/blind test splits 
#' (i.e. iterations of outer loop)
#' @param numCores (integer) number of CPU cores for parallel processing
#' @param JavaMemory (integer) memory in (Gb) used for each fold of CV
#' @param featSelCutoff (integer) cutoff for inner-fold CV to call 
#' feature-selected in a given split
#' @param keepAllData (logical) if TRUE keeps all intermediate files, even
#' those not needed for assessing the predictor. Use very cautiously as for
#' some designs, each split can result in using 1Gb of data.
#' @param startAt (integer) which of the splits to start at (e.g. if the
#' job aborted part-way through)
#' @param preFilter (logical) if TRUE uses lasso to prefilter dataList within 
#' cross-validation loop. Only variables that pass lasso get included. The
#' current option is not recommended for pathway-level features as most genes
#' will be eliminated by lasso. Future variations may allow other prefiltering
#' options that are more lenient.
#' @param preFilterGroups (char) vector with subset of names(dataList)
#' to which prefiltering needs to be limited. Allows users to indicate
#' which data layers should be prefiltered using regression and which
#' are to be omitted from this process. Prefiltering uses regression, which
#' omits records with missing values. Structured missingness can result in
#' empty dataframes if missing values are removed from these, which in turn
#' can crash the predictor. To impute missing data, see the 'impute' and 
#' 'imputeGroups' parameters. 
#' @param impute (logical) if TRUE applies imputation by median within CV
#' @param imputeGroups (char) If impute set to TRUE, indicate which groups you 
#' want imputed. 
#' @param debugMode (logical) when TRUE runs jobs in serial instead of parallel and 
#' prints verbose messages. Also prints system Java calls and prints all standard out
#' and error output associated with these calls.
#' @param logging (char) level of detail with which messages are printed. 
#' Options are: 1) none: turn off all messages; 2) all: greatest level of 
#' detail (recommended for advanced users, or for debugging); 3) default: 
#' print key details (useful setting for most users)
#' @import glmnet
#' @importFrom stats median na.omit coef
#' @importFrom utils read.delim write.table
#' @importFrom methods is
#' @return (list) 
#' "inputNets": data.frame of all input network names. Columns are "NetType"
#' (group) and "NetName" (network name).
#' "Split<i>" is the data for train/test split i 
#' (i.e. one per train/test split).
#' Each "SplitX" entry contains in turn a list of results for that split. 
#' Key-value pairs are:
#' 1) predictions: real and predicted labels for test patients
#' 2) accuracy: percent accuracy of predictions
#' 3) featureScores: list of length g, where g is number of patient classes.
#' scores for all features following feature selection, for corresponding 
#' class.
#' 4) featureSelected: list of length g (num patient classes). List of 
#' selected features for corresponding patient class, for that train/test 
#' split. Side effect of generating predictor-related data in <outDir>.
#' @export
#' @importFrom methods is
#' @import MultiAssayExperiment
#' @examples
#'
#' library(curatedTCGAData)
#' library(MultiAssayExperiment)
#' curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)
#' 
#' # fetch mrna, mutation data
#' brca <- curatedTCGAData("BRCA",c("mRNAArray"),FALSE)
#' 
#' # get subtype info
#' pID <- colData(brca)$patientID
#' pam50 <- colData(brca)$PAM50.mRNA
#' staget <- colData(brca)$pathology_T_stage
#' st2 <- rep(NA,length(staget))
#' st2[which(staget %in% c("t1","t1a","t1b","t1c"))] <- 1
#' st2[which(staget %in% c("t2","t2a","t2b"))] <- 2
#' st2[which(staget %in% c("t3","t3a"))] <- 3
#' st2[which(staget %in% c("t4","t4b","t4d"))] <- 4
#' pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"                         
#' pam50[which(pam50 %in% "Luminal A")] <- "LumA"
#' colData(brca)$ID <- pID
#' colData(brca)$STAGE <- st2                                                 
#' colData(brca)$STATUS <- pam50
#' 
#' # keep only tumour samples
#' idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
#' cat(sprintf("excluding %i samples\n", length(idx)))
#'                                                                            
#' tokeep <- setdiff(pID, pID[idx])
#' brca <- brca[,tokeep,]
#' 
#' pathList <- readPathways(fetchPathwayDefinitions())
#' brca <- brca[,,1] # keep only clinical and mRNA data
#' 
#' # remove duplicate arrays
#' smp <- sampleMap(brca)
#' samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
#' notdup <- samps[which(!duplicated(samps$primary)),"colname"]
#' brca[[1]] <- brca[[1]][,notdup]
#' 
#' groupList <- list()
#' groupList[["BRCA_mRNAArray-20160128"]] <- pathList[seq_len(3)]
#' groupList[["clinical"]] <- list(
#'	age="patient.age_at_initial_pathologic_diagnosis",
#'  stage="STAGE")
#' makeNets <- function(dataList, groupList, netDir,...) {
#'     netList <- c()
#'     # make RNA nets: group by pathway
#'     if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) {
#'     netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
#'                 rownames(dataList[["BRCA_mRNAArray-20160128"]]),
#'                 groupList[["BRCA_mRNAArray-20160128"]],
#'                 netDir,verbose=FALSE,
#'                 writeProfiles=TRUE,...)
#'     netList <- unlist(netList)
#'     cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
#'     }
#' 
#'     # make clinical nets,one net for each variable
#'     netList2 <- c()
#'     if (!is.null(groupList[["clinical"]])) {
#'     netList2 <- makePSN_NamedMatrix(dataList$clinical,
#'         rownames(dataList$clinical),
#'         groupList[["clinical"]],netDir,
#'         simMetric="custom",customFunc=normDiff, # custom function
#'         writeProfiles=FALSE,
#'         sparsify=TRUE,verbose=TRUE,...)
#'     }
#'     netList2 <- unlist(netList2)
#'     cat(sprintf("Made %i clinical nets\n", length(netList2)))
#'     netList <- c(netList,netList2)
#'     cat(sprintf("Total of %i nets\n", length(netList)))
#'     return(netList)
#' }
#' 
#' # takes 10 minutes to run
#' #out <- buildPredictor(dataList=brca,groupList=groupList,
#' #   makeNetFunc=makeNets, ### custom network creation function
#' #   outDir=paste(tempdir(),"pred_output",sep=.Platform$file.sep), ## absolute path
#' #   numCores=16L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)
buildPredictor <- function(dataList,groupList,outDir=tempdir(),makeNetFunc,
	featScoreMax=10L,trainProp=0.8,numSplits=10L,numCores,JavaMemory=4L,
	featSelCutoff=9L,keepAllData=FALSE,startAt=1L, preFilter=FALSE,
	impute=FALSE,preFilterGroups=NULL, imputeGroups=NULL,logging="default",
	debugMode=FALSE) { 
verbose_default <- TRUE
verbose_runQuery <- FALSE	  # messages when running individual queries
verbose_compileNets <- FALSE  # message when compiling PSN into database
verbose_runFS <- TRUE		  # runFeatureSelection() 
verbose_predict <- FALSE
verbose_compileFS <- FALSE
verbose_makeFeatures <- FALSE

if (logging == "all") {
	verbose_runQuery <- TRUE
	verbose_compileNets <- TRUE 
	verbose_compileFS <- TRUE
	verbose_makeFeatures <- TRUE
} else if (logging=="none") {
	verbose_runFS<-FALSE
	verbose_default <- FALSE
	verbose_predict <- FALSE
}

# Check input
if (missing(dataList)) stop("dataList must be supplied.\n")
if (missing(groupList)) stop("groupList must be supplied.\n")
if (length(groupList)<1) stop("groupList must be of length 1+\n")
tmp <- unlist(lapply(groupList,class))
not_list <- sum(tmp == "list")<length(tmp)
nm1 <-setdiff(names(groupList),"clinical") 
if (!is(dataList,"MultiAssayExperiment"))
	stop("dataList must be a MultiAssayExperiment")
names_nomatch <- any(!nm1 %in% names(dataList))
if (!is(groupList,"list") || not_list || names_nomatch ) {
	msg <- c("groupList must be a list of lists.",
	" Names must match those in dataList, and each entry should be a list",
  " of networks for this group.")
	stop(paste(msg,sep=""))
}

if (outDir != normalizePath(outDir)) {
	browser()
	stop("outDir should be an absolute path, not relative.")
}

if (trainProp <= 0 | trainProp >= 1) 
		stop("trainProp must be greater than 0 and less than 1")
if (startAt > numSplits) stop("startAt should be between 1 and numSplits")

megaDir <- outDir
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

# set aside for testing within each split
pheno_all <- colData(dataList)
pheno_all <- as.data.frame(pheno_all)

message("Predictor started at:")
message(Sys.time())

# run featsel once per subtype
subtypes <- unique(pheno_all$STATUS)

# convert to list structure 
exprs <- experiments(dataList)
datList2 <- list()
for (k in seq_len(length(exprs))) {
	tmp <- exprs[[k]]
	df <- sampleMap(dataList)[which(sampleMap(dataList)$assay==names(exprs)[k]),]
	colnames(tmp) <- df$primary[match(df$colname,colnames(tmp))]
	tmp <- as.matrix(assays(tmp)[[1]]) # convert to matrix
	datList2[[names(exprs)[k]]]<- tmp	
}
if ("clinical" %in% names(groupList)) {
	tmp <- colData(dataList)
	vars <- unique(unlist(groupList[["clinical"]]))
	datList2[["clinical"]] <- t(as.matrix(tmp[,vars,drop=FALSE]))
}

dataList <- datList2; rm(datList2);

if (verbose_default){
	message(sprintf("-------------------------------"))
	message(sprintf("# patients = %i", nrow(pheno_all)))
	message(sprintf("# classes = %i { %s }", length(subtypes),
		paste(subtypes,collapse=",")))
	message("Sample breakdown by class")
	message(table(pheno_all$STATUS))
	message(sprintf("%i train/test splits",numSplits))
	message(sprintf("Feature selection cutoff = %i of %i",
		featSelCutoff,featScoreMax))
	message(sprintf("Datapoints:"))
	for (nm in names(dataList)) {
		message(sprintf("\t%s: %i units", nm, nrow(dataList[[nm]])))
	}
}


outList <- list()

# create master list of possible networks
tmp <- list()
for (nm in names(groupList)) {
	curNames <- names(groupList[[nm]])
	tmp[[nm]] <- cbind(rep(nm,length(curNames)),curNames)
}
tmp <- do.call("rbind",tmp)
if (length(nm) < 2) tmp <- as.matrix(tmp)
colnames(tmp) <- c("NetType","NetName")
outList[["inputNets"]] <- tmp

if (verbose_default) {
	message("\n\nCustom function to generate input nets:")
	print(makeNetFunc)
	message(sprintf("-------------------------------\n"))
}

for (rngNum in startAt:numSplits) {
	curList <- list()

	if (verbose_default) {
	message(sprintf("-------------------------------"))
	message(sprintf("Train/test split # %i", rngNum))
	message(sprintf("-------------------------------"))
	}
	outDir <- paste(megaDir,sprintf("rng%i",rngNum),sep=.Platform$file.sep)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
		verbose=verbose_default)
	pheno <- pheno_all[which(pheno_all$TT_STATUS %in% "TRAIN"),]

	dats_train <- lapply(dataList, function(x) 
		x[,which(colnames(x) %in% pheno$ID)])

	if (impute) {
	if (verbose_default) message("**** IMPUTING ****")
	if (is.null(imputeGroups)) imputeGroups <- names(dats_train)
	if (!any(imputeGroups %in% names(dats_train))) 
		stop("imputeGroups must match names in dataList")

	nmset <- names(dats_train)
	dats_train <- lapply(names(dats_train), function(nm) {
		x <- dats_train[[nm]]
		print(class(x))
		if (nm %in% imputeGroups) {
			missidx <- which(rowSums(is.na(x))>0) 
			for (i in missidx) {
				na_idx <- which(is.na(x[i,]))
				x[i,na_idx] <- median(x[i,],na.rm=TRUE) 
			}
		} 
		x
	})
	names(dats_train) <- nmset
	}

	# prefilter with lasso
	if (preFilter) {
	if (is.null(preFilterGroups)) preFilterGroups <- names(dats_train)
	if (!any(preFilterGroups %in% names(dats_train))) {
		stop("preFilterGroups must match names in dataList")
	}
	
	message("Prefiltering enabled")
	for (nm in preFilterGroups) {
		message(sprintf("%s: %i variables",nm,nrow(dats_train[[nm]])))
		if (nrow(dats_train[[nm]])<2)  # only has one var, take it.
			vars <- rownames(dats_train[[nm]])
		else { 
			newx <- na.omit(dats_train[[nm]])
			tmp <- pheno[which(pheno$ID %in% colnames(newx)),]
			tryCatch( {
			fit <- cv.glmnet(x=t(newx),
					y=factor(tmp$STATUS), family="binomial", alpha=1) # lasso
			}, error=function(ex) {
				print(ex)
				message("*** You may need to set impute=TRUE for prefiltering ***")
			},finally={
			})
			wt <- abs(coef(fit,s="lambda.min")[,1])
			vars <- setdiff(names(wt)[which(wt>.Machine$double.eps)],
				"(Intercept)")
			}
		if (length(vars)>0) {
			tmp <- dats_train[[nm]]
			tmp <- tmp[which(rownames(tmp) %in% vars),,drop=FALSE]
			dats_train[[nm]] <- tmp
		} else {
			# leave dats_train as is, make a single net
		} 
		message(sprintf("rngNum %i: %s: %s pruned",rngNum,nm,length(vars)))
		}
	}
	
	if (verbose_default) {
		message("# values per feature (training)")
		for (nm in names(dats_train)) {
			message(sprintf("\tGroup %s: %i values", 
				nm,nrow(dats_train[[nm]])))
		}
	}

	netDir <- paste(outDir,"tmp",sep=.Platform$file.sep)
	dir.create(netDir)
message("about to setup featuredb")
	pheno_id <- setupFeatureDB(pheno,netDir)
message("done setting up feature db")

	if (verbose_default) message("** Creating features")
	createPSN_MultiData(dataList=dats_train,groupList=groupList,
			pheno=pheno_id,
			netDir=netDir,customFunc=makeNetFunc,numCores=numCores,
			verbose=verbose_makeFeatures)
	if (verbose_default) message("** Compiling features")
	dbDir <- compileFeatures(netDir,outDir, numCores=numCores, 
			verbose=verbose_compileFS, debugMode=debugMode)
	if (verbose_default) message("\n** Running feature selection")

	curList[["featureScores"]] <- list()

	for (g in subtypes) {
	    pDir <- paste(outDir,g,sep=.Platform$file.sep)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE);
			dir.create(pDir)
			if (verbose_default) message(sprintf("\tClass: %s",g))
			pheno_subtype <- pheno
			pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
			trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% g)]
			if (verbose_default) {
				print(table(pheno_subtype$STATUS,useNA="always"))
}
		
			# Cross validation
			resDir <- paste(pDir,"GM_results",sep=.Platform$file.sep)
			message(sprintf("\tScoring features"))
			runFeatureSelection(trainPred, 
				outDir=resDir, dbPath=dbDir$dbDir, 
				nrow(pheno_subtype),verbose=verbose_runFS, 
				numCores=numCores, verbose_runQuery=TRUE, # verbose_runQuery,
				featScoreMax=featScoreMax,JavaMemory=JavaMemory,
				debugMode=debugMode)
	
	  	# Compute network score
			nrank <- dir(path=resDir,pattern="NRANK$")
			if (verbose_default) message("\tCompiling feature scores")
			pTally <- compileFeatureScores(paste(resDir,nrank,
					sep=.Platform$file.sep),
				verbose=verbose_compileFS)
			tallyFile <- paste(resDir,
				sprintf("%s_pathway_CV_score.txt",g),
				sep=.Platform$file.sep)
			write.table(pTally,file=tallyFile,sep="\t",
				col.names=TRUE,row.names=FALSE,
				quote=FALSE)
			curList[["featureScores"]][[g]] <- pTally
		if (verbose_default) message("")
	}
	
	## Class prediction for this split
	if (verbose_default) message("\n** Predicting labels for test")
	pheno <- pheno_all
	predRes <- list()

	curList[["featureSelected"]] <- list()
	for (g in subtypes) {
		if (verbose_default) message(sprintf("%s",g))
		pDir <- paste(outDir,g,sep=.Platform$file.sep)
		pTally <- read.delim(
			paste(pDir,"GM_results",
				sprintf("%s_pathway_CV_score.txt",g),
				sep=.Platform$file.sep),
			sep="\t",header=TRUE,as.is=TRUE)
		idx <- which(pTally[,2]>=featSelCutoff)

		pTally <- pTally[idx,1]
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont.txt","",pTally)

		curList[["featureSelected"]][[g]] <- pTally

		if (verbose_default)
			message(sprintf("\t%i feature(s) selected",length(pTally)))
		netDir <- paste(pDir,"networks",sep=.Platform$file.sep)

		dats_tmp <- list()
		for (nm in names(dataList)) {
			passed <- rownames(dats_train[[nm]])
			tmp <- dataList[[nm]]
			# only variables passing prefiltering should be used to make PSN
			dats_tmp[[nm]] <- tmp[which(rownames(tmp) %in% passed),,
				drop=FALSE] 
		}		

		# ------
		# Impute test samples if flag set
		# impute
		if (impute) {
			train_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TRAIN")]
			test_samp <- pheno_all$ID[which(pheno_all$TT_STATUS %in% "TEST")]
			nmSet <- names(dats_tmp)
			dats_tmp <- lapply(names(dats_tmp), function(nm) {
				x <- dats_tmp[[nm]]
				if (nm %in% imputeGroups) {
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
				}
				x
			})
			names(dats_tmp) <- nmSet
			#alldat_tmp <- do.call("rbind",dats_tmp)
			}

		if (verbose_default) message(sprintf("\tCreate & compile features",g))
		if (length(pTally)>=1) {
		netDir <- paste(pDir,"tmp",sep=.Platform$file.sep)
		dir.create(netDir)
		pheno_id <- setupFeatureDB(pheno,netDir)
		createPSN_MultiData(dataList=dats_tmp,groupList=groupList,
			pheno=pheno_id,
			netDir=netDir,customFunc=makeNetFunc,numCores=numCores,
			filterSet=pTally,verbose=verbose_default)
		dbDir <- compileFeatures(netDir,outDir=pDir,numCores=numCores,
			verbose=verbose_compileNets,debugMode=debugMode)

		# run query for this class
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
		qFile <- paste(pDir,sprintf("%s_query",g),sep=.Platform$file.sep)
		writeQueryFile(qSamps,"all",nrow(pheno),qFile)
		if (verbose_default) message(sprintf("\t** %s: Compute similarity",g))
		resFile <- runQuery(dbDir$dbDir,qFile,resDir=pDir,
			JavaMemory=JavaMemory, numCores=numCores,
			verbose=verbose_runQuery,debugMode=debugMode)
		predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK",resFile),pheno,g)
		} else {
			predRes[[g]] <- NA
		}
	}
	if (verbose_default) message("")
	
	if (sum(is.na(predRes))>0 & verbose_default) {
		str <- sprintf("RNG %i : One or more classes have no selected features.",
				rngNum)
		str <- sprintf("%s Not classifying.",str)
		message(str)
	} else {
		if (verbose_default) message("** Predict labels")
		predClass <- predictPatientLabels(predRes,
			verbose=verbose_predict)
		out <- merge(x=pheno_all,y=predClass,by="ID")
		outFile <- paste(outDir,"predictionResults.txt",
			sep=.Platform$file.sep)
		acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
		if (verbose_default)
			message(sprintf("Split %i: ACCURACY (N=%i test) = %2.1f%%",
			rngNum, nrow(out), acc*100))

		curList[["predictions"]] <- out
		curList[["accuracy"]] <- acc
	}
        
	if (!keepAllData) {
		unlink(outDir, recursive=TRUE)
	}# endif !keepAllData

	if (verbose_default) {
		message("\n----------------------------------------")
	}

	outList[[sprintf("Split%i",rngNum)]] <- curList
	}
	message("Predictor completed at:")
	message(Sys.time())

	return(outList)
}

