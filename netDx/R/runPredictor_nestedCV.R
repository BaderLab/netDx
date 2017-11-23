#' Run nested cross-validation on data
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
#' See examples/NestedCV_MultiData.Rmd for a full working example.
#' @param pheno (data.frame) sample metadata, must have ID and STATUS columns
#' @param dataList (list) keys are datatypes; values contain patient data
#' for the corresponding datatype. e.g. dataList[["rna"]] contains expression
#' matrix. Rows are units (e.g. genes, individual clinical variables) and
#' columns are patients
#' @param groupList (list) keys are datatypes and values are lists indicating
#' how units for those datatypes are to be grouped. Keys must match those
#' in dataList. Each entry of groupList[[k]] will generate a new PSN.
#'  e.g. groupList[["rna"]] could be a list of pathway definitions.
#' So keys(groupList[["rna"]]) would have pathway names, generating one PSN
#' per pathways, and values(groupList[["rna"]]) would be genes that would be
#' grouped for the corresponding pathwayList.
#' @param makeNetFunc (function) user-defined function for creating the set
#' of input PSN provided to netDx. See createPSN_MultiData()::customFunc.
#' @param outDir (char) directory where results will be stored. If this
#' directory exists, its contents will be overwritten
#' @param trainProp (numeric 0 to 1) Percent samples to use for training
#' @param nFoldCV (integer) number of CV folds in inner loop
#' @param numSplits (integer) number of train/blind test splits (i.e. iterations
#' of outer loop)
#' @param numCores (integer) number of CPU cores for parallel processing
#' @param CVmemory (integer) memory in (Gb) used for each fold of CV
#' @param CVcutoff (integer) cutoff for inner-fold CV to call feature-selected
#' in a given split
#' @param keepAllData (logical) if TRUE keeps all intermediate files, even
#' those not needed for assessing the predictor. Use very cautiously as for
#' some designs, each split can result in using 1Gb of data.
#' @examples see examples/NestedCV_MultiData.Rmd for example use.
#' @export

runPredictor_nestedCV <- function(pheno,dataList,groupList,outDir,makeNetFunc,
	typeList=NULL,nFoldCV=10L,trainProp=0.8,numSplits=10L,
	numCores,CVmemory=4L,CVcutoff=9L,cliqueReps=2500L,cliquePthresh=0.07,keepAllData=FALSE) {

### tests# pheno$ID and $status must exist
if (missing(dataList)) stop("dataList must be supplied.\n")
if (missing(groupList)) stop("groupList must be supplied.\n")
if (missing(typeList)) stop("typeList must be supplied.\n")
if (trainProp <= 0 | trainProp >= 1)
		stop("trainProp must be greater than 0 and less than 1")

megaDir <- outDir
if (file.exists(megaDir)) unlink(megaDir,recursive=TRUE)
dir.create(megaDir)

# set aside for testing within each split
pheno_all <- pheno;

#track clique status
stat_mat <- list()

num_preds <- 0

#If typeList not given, make it assuming all continuous data types
if(is.null(typeList)){
	typeList <- list()
	for(cur_type in names(dataList)){
		typeList[[cur_type]] <- 'continuous'
	}
}

logFile <- sprintf("%s/log.txt",megaDir)
sink(logFile,split=TRUE)
cat("Predictor started at:\n")
print(Sys.time())
tryCatch({

# run featsel once per subtype
subtypes <- unique(pheno$STATUS)

cat(sprintf("-------------------------------\n"))
cat(sprintf("# patients = %i\n", nrow(pheno)))
cat(sprintf("# classes = %i { %s }\n", length(subtypes),
	paste(subtypes,collapse=",")))
cat("Sample breakdown by class\n")
print(table(pheno$STATUS))
cat(sprintf("Nested CV design = %i CV x %i splits\n", nFoldCV, numSplits))
cat(sprintf("Datapoints:\n"))
for (nm in names(dataList)) {
	cat(sprintf("\t%s: %i units\n", nm, nrow(dataList[[nm]])))
}

# create master list of possible networks
cat("# input nets provided:\n")
netFile <- sprintf("%s/inputNets.txt", megaDir)
cat("NetType\tNetName\n",file=netFile)
for (nm in names(groupList)) {
	curNames <- names(groupList[[nm]])
	for (nm2 in curNames) {
		cat(sprintf("%s\t%s\n",nm,nm2),file=netFile,append=TRUE)
	}
}

cat("\n\nCustom function to generate input nets:\n")
print(makeNetFunc)
cat(sprintf("-------------------------------\n\n"))

for (rngNum in 1:numSplits) {
	stat_mat[[rngNum]] <- list(rng = rngNum,
															binNets_preCFilter = '',
															binNets_postCFilter = '',
															nonBinNets = '',
															patients_preCFilter = '',
															patients_postCFilter = ''
	)
	cat(sprintf("-------------------------------\n"))
	cat(sprintf("RNG seed = %i\n", rngNum))
	cat(sprintf("-------------------------------\n"))
	outDir <- sprintf("%s/rng%i",megaDir,rngNum)
	dir.create(outDir)

	pheno_all$TT_STATUS <- splitTestTrain(pheno_all,pctT=trainProp,
											  setSeed=rngNum*5)
	pheno <- subset(pheno_all, TT_STATUS %in% "TRAIN")

	#Extract training samples from data
	temp_list <- list()
	for(cur_dat in names(dataList)){
		if(class(dataList[[cur_dat]])[1] == 'GRanges'){
			dats_train_temp <- dataList[[cur_dat]][mcols(dataList[[cur_dat]])$ID %in% pheno$ID]
			temp_list[[cur_dat]] <- dats_train_temp
		}else{
			dats_train_temp <- dataList[[cur_dat]][,which(colnames(dataList[[cur_dat]]) %in% pheno$ID)]
			temp_list[[cur_dat]] <- dats_train_temp
			}
	}
	dats_train <- temp_list

	netDir <- sprintf("%s/networks",outDir)
	netList <- createPSN_MultiData(dataList=dats_train,groupList=groupList,
			netDir=netDir,customFunc=makeNetFunc,numCores=numCores,typeList=typeList)

	#take netList and update net names for binary data types according to typeList
	bin_netList <- c()
	nonBin_netList <- c()
	for(cur_name in names(typeList)){
		if(typeList[[cur_name]] == 'binary'){
			netList[[cur_name]] <- updateBinaryNetNames(netDir,netList[[cur_name]])
			bin_netList <- c(bin_netList, netList[[cur_name]])
		}else{
			nonBin_netList <- c(nonBin_netList, netList[[cur_name]])
		}
	}

	stat_mat[[rngNum]][['binNets_preCFilter']] <- length(bin_netList)
	stat_mat[[rngNum]][['nonBinNets']] <- length(nonBin_netList)

	if(length(bin_netList) >= 1){
		tmp <- updateBinaryNetlist(netDir,bin_netList, pheno)
		p_train_og <- tmp[[1]]
		pheno_train_og <- tmp[[2]]
	}else{
		pheno_train_og <- pheno
		}

	for (g in subtypes) {
	    pDir <- sprintf("%s/%s",outDir,g)
	    if (file.exists(pDir)) unlink(pDir,recursive=TRUE);dir.create(pDir)

			#Run clique-filtering only for binary nets
			if(length(bin_netList) >= 1){
				temp_netdir <- paste(netDir,"_bin", sep = '')
				if (!file.exists(temp_netdir)){
					dir.create(temp_netdir)
					for(cur_file in bin_netList){
						file.rename(from = sprintf('%s/%s',netDir, cur_file),
												to = sprintf('%s/%s',temp_netdir, cur_file))
					}
				}else{
					for(cur_file in bin_netList){
						#If we have done one subtype, now we must delete the clique networks
						#from the netDir
						if (file.exists(sprintf('%s/%s',netDir, cur_file))){
							unlink(sprintf('%s/%s',netDir, cur_file))
						}
					}
				}

				cat("Running clique-filtering\n")
				netInfo <- cliqueFilterNets(temp_netdir,pheno_train_og,outDir,
					predClass=g,numReps=cliqueReps,numCores=numCores)
				pvals   <- as.numeric(netInfo[,"pctl"])

				netInfo <- netInfo[which(pvals < cliquePthresh),]
				print(nrow(netInfo))
				p_train  	<- p_train_og[,
					which(colnames(p_train_og) %in% rownames(netInfo))]

				stat_mat[[rngNum]][['binNets_postCFilter']] <- nrow(netInfo)

				# update nets after clique-filtering
				cat("Clique filtered\n")
				trainNetDir <- sprintf("%s/%s_networksCliqueFilt",outDir,g)
				tmp <- updateNets(p_train, pheno_train_og,
								oldNetDir=temp_netdir, newNetDir=trainNetDir)
				p	<- tmp[[1]]
				pheno_clique	<- tmp[[2]]

				#COPY NETWROKS FROM CLIQUE OVER TO NETWORKS FOLDER THEN USE THAT AS YOUR
				#NETDIR
				for(cur_file in rownames(netInfo)){
					if(file.exists(sprintf('%s/%s',trainNetDir, cur_file))){
						file.rename(from = sprintf('%s/%s',trainNetDir, cur_file),
												to = sprintf('%s/%s',netDir, cur_file))
						}
					}

				unlink(trainNetDir, recursive=TRUE)

				stat_mat[[rngNum]][['patients_preCFilter']] <- nrow(pheno_train_og)
				stat_mat[[rngNum]][['patients_postCFilter']] <- nrow(pheno_clique)

			}else{
				#If there are no binary nets, skip clique filtering
				pheno_train <- pheno_train_og
			}

			# Clqiue filtering is only run on binary netwroks, thus some genetic data
			# patients are excluded from the resulting pheno data. Append them back in
			# only if clique filtering was run
			if((length(bin_netList) >= 1) & (length(nonBin_netList) >= 1)){
				for(cur_datatype in names(typeList)){
					if(typeList[[cur_datatype]] != 'binary'){
						cur_phenoIDs <- colnames(dataList[[cur_datatype]])
						temp_pheno <- pheno_all[which(cur_phenoIDs %in% pheno_all$ID),]
						temp_pheno <- temp_pheno[which(temp_pheno$TT_STATUS %in% 'TRAIN'),]
						pheno_train <- rbind(pheno_clique, temp_pheno)
					}
				}
			}else if((length(bin_netList) >= 1) & (length(nonBin_netList) == 0)){
				pheno_train <- pheno_clique
				}

			#Generate GeneMANIA database
			dbDir	<- GM_createDB(netDir, pheno_train$ID, outDir,numCores=numCores)
			cat(sprintf("\n******\nClass: %s\n",g))
			pheno_subtype <- pheno_train
			pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
			trainPred <- pheno_subtype$ID[which(pheno_subtype$STATUS %in% g)]
			print(table(pheno_subtype$STATUS,useNA="always"))

			# Cross validation
			resDir <- sprintf("%s/GM_results",pDir)
			GM_runCV_featureSet(trainPred,
				outDir=resDir, GM_db=dbDir$dbDir,
				nrow(pheno_subtype),verbose=T, numCores=numCores,
				nFold=nFoldCV,GMmemory=CVmemory)

	  	# Compute network score
			nrank <- dir(path=resDir,pattern="NRANK$")
			pTally		<- GM_networkTally(paste(resDir,nrank,sep="/"))
			tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
			write.table(pTally,file=tallyFile,sep="\t",col=T,row=F,quote=F)
	}

	#clean up unused variables
	rm(pheno_train_og)
	rm(pheno_clique)

	## Class prediction for this split
	pheno <- pheno_all
	predRes <- list()
	no_pred <- 0
	for (g in subtypes) {
		pDir <- sprintf("%s/%s",outDir,g)
		pTally <- read.delim(
			sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir,g),
			sep="\t",h=T,as.is=T)
		pTally <- pTally[which(pTally[,2]>=CVcutoff),1]
		pTally <- sub(".profile","",pTally)
		pTally <- sub("_cont","",pTally)
		pTally <- sub("binary_","",pTally)
		if(identical(pTally, character(0))){
			print("NO PTALLY!!!!!!!!!")
			no_pred <- 1
			next;
		}
		cat(sprintf("%s: %i networks\n",g,length(pTally)))
		netDir <- sprintf("%s/networks",pDir)

		netList<- createPSN_MultiData(dataList=dataList,groupList=groupList,
			netDir=netDir,
			customFunc=makeNetFunc,numCores=numCores, typeList=typeList,
			filterSet=pTally)

		dbDir <- GM_createDB(netDir,pheno$ID,pDir,numCores=numCores)

		# run query for this class
		qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
		qFile <- sprintf("%s/%s_query",pDir,g)
		GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
		resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)
		predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
		#doesn't work if we do not reset pheno
		pheno <- pheno_all
	}
 	if(no_pred != 1){
		predClass <- GM_OneVAll_getClass(predRes)
		out <- merge(x=pheno_all,y=predClass,by="ID")
		outFile <- sprintf("%s/predictionResults.txt",outDir)
		write.table(out,file=outFile,sep="\t",col=T,row=F,quote=F)

		acc <- sum(out$STATUS==out$PRED_CLASS)/nrow(out)
		cat(sprintf("Accuracy on %i blind test subjects = %2.1f%%\n",
			nrow(out), acc*100))
		num_preds <- num_preds + 1
	}else{
		cat(sprintf("no prediction made for rng %s, no feat sel nets", rngNum))
	}

	if (!keepAllData) {
    system(sprintf("rm -r %s/dataset %s/tmp %s/networks %s/networks_bin",
        outDir,outDir,outDir, outDir))
	for (g in subtypes) {
    system(sprintf("rm -r %s/%s/dataset %s/%s/networks",
        outDir,g,outDir,g))
	}
	}}
	num_pred_total <- (num_preds/numSplits)*100
	cat(sprintf("%s percent of splits had predictions generated \n",num_pred_total))
	stat_df <- do.call(rbind, lapply(stat_mat, data.frame, stringsAsFactors=FALSE))
	write.csv(stat_df, file = sprintf("%s/networkInfo.csv",megaDir), quote = FALSE)
}, error=function(ex){
	print(ex)
}, finally={
	cat("Predictor completed at:\n")
	print(Sys.time())
	sink(NULL)
})

}
