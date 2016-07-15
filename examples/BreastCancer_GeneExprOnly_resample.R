## ----set-params,cache=FALSE,eval=TRUE------------------------------------
rm(list=ls())

# Change this to a local directory where you have write permission
outDir <- "~/tmp/TCGA_BRCA_geneXpr_resample" 

numCores 	<- 8L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
cutoff		<- 9L  	# score cutoff for feature-selected networks
TRAIN_PROP <- 0.67 	# fraction of samples to use for training
numResamples <- 3	# number of data resamplings

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

## ----load-packages, cache=FALSE,eval=TRUE--------------------------------
require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

sink("BreastCancer_GeneExprOnly.log",split=TRUE)
tryCatch({

## ----split-test-train, cache=FALSE,eval=TRUE-----------------------------

subtypes<- c("LumA")
pheno$STATUS[which(!pheno$STATUS %in% subtypes)] <- "other"
subtypes <- c(subtypes,"other") # add residual

pheno$TT_STATUS <- splitTestTrain(pheno,
    pctT = TRAIN_PROP,setSeed = 42,predClass = "LumA" )

## ----clean-pheno-xpr, cache=FALSE,eval=TRUE------------------------------
pheno_FULL	<- pheno
xpr_FULL 	<- xpr
###cnv_FULL	<- cnv_GR
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
xpr			<- xpr[,which(colnames(xpr)%in% pheno$ID)]
###cnv_GR		<- cnv_GR[which(cnv_GR$ID %in% pheno$ID)]

## ----read-pathways, cache=FALSE,eval=TRUE--------------------------------
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
    path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
head(pathwayList)

## ----make-xpr-psn, eval=TRUE---------------------------------------------
###profDir <- sprintf("%s/profiles",outDir)
###netDir <- sprintf("%s/networks",outDir)

###netList <- makePSN_NamedMatrix(xpr, rownames(xpr), 
###        pathwayList,profDir,verbose=FALSE,
###        numCores=numCores,writeProfiles=TRUE)
###netList <- unlist(netList)
###head(netList)

## ----make-cnv-psn, cache=FALSE,eval=TRUE---------------------------------
data(genes)
gene_GR     <- GRanges(genes$chrom,
   IRanges(genes$txStart,genes$txEnd),
   name=genes$name2)

## ----feature-selection,cache=FALSE,eval=TRUE-----------------------------
## repeat process for each class
for (g in subtypes) {
    pDir <- sprintf("%s/%s",outDir,g)
    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
	dir.create(pDir)

	cat(sprintf("\n******\nSubtype %s\n",g))
	pheno_subtype <- pheno

	### run three resamplings, with 70% training in each.
	### run 10-fold cross validation using the training population.
	### finally add up the score across the three resamplings
	### this score will be used to assess performance on the test samples
	### put aside at the very beginning.
	for (repNum in 1:numResamples) {
		cat(sprintf("Resampling %i ***********\n",repNum))
		pheno_subtype <- pheno
		pDir <- sprintf("%s/%s/part%i", outDir,g,repNum)
    	if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
		dir.create(pDir)

		### TODO Change this to use splitter that uses one patient for test
		### exactly once.
		pheno_subtype$INNER_TT_STATUS <- splitTestTrain(pheno_subtype,
	    	pctT = TRAIN_PROP,predClass = g )
		pheno_subtype <- subset(pheno_subtype, INNER_TT_STATUS%in% "TRAIN")
		cat(sprintf("Resampling %i: %i training ***\n" ,repNum, 
					nrow(pheno_subtype)))
		
		## label patients not in the current class as a residual
		pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)]<- "nonpred"
		print(table(pheno_subtype$STATUS,useNA="always"))

		# prepare nets
		profDir <- sprintf("%s/profiles",pDir)
		netDir <- sprintf("%s/networks",pDir)
		xpr_tmp	<- xpr[,which(colnames(xpr)%in% pheno_subtype$ID)]
		print(dim(xpr_tmp))
		netList <- makePSN_NamedMatrix(xpr_tmp, rownames(xpr_tmp), 
		        pathwayList,profDir,verbose=FALSE,
		        numCores=numCores,writeProfiles=TRUE)
		netList <- unlist(netList)
		head(netList)

		## prune nets to limit to patients in this set
		dbDir	<- GM_createDB(profDir, pheno_subtype$ID, outDir,
							numCores=numCores)
		resDir	<- sprintf("%s/GM_results",pDir)
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

	## finally tally across resamplings and write final result file
	pathScore <- list()
	for (currRep in 1:numResamples) {
		f <- sprintf("%s/%s/part%i/GM_results/%s_pathway_CV_score.txt",
			outDir,g,currRep,g)
		dat <- read.delim(f,sep="\t",header=TRUE,as.is=TRUE)
		for (k in 1:nrow(dat)) {
			p <- dat[k,1]
			if (p %in% names(pathScore)) 
					pathScore[[p]] <-pathScore[[k]]+dat[k,2]
			else pathScore[[p]] <- dat[k,2]
		}
	}
	df <- data.frame(PATHWAY_NAME=names(pathScore), SCORE=unlist(pathScore))
	write.table(df,file=sprintf("%s/%s/%s_pathwayScore.txt", outDir,
			g,g),sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}
cat("DONE feature selection!\n")

## ----class-prediction, eval=TRUE-----------------------------------------
# now create GM databases for each class
# should contain train + test patients
# and be limited to nets that pass feature selection
pheno <- pheno_FULL
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
	profDir <- sprintf("%s/profiles",pDir)

	# prepare nets for new db
	tmp <- makePSN_NamedMatrix(xpr_FULL,rownames(xpr),
		pathwayList[which(names(pathwayList)%in% pTally)],
		profDir,verbose=F,numCores=numCores,writeProfiles=TRUE)
	###tmp <- makePSN_RangeSets(cnv_FULL,
	###	path_GRList[which(names(path_GRList)%in% pTally)],
	###		profDir,verbose=FALSE)
	# create db
	dbDir <- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

	# query of all training samples for this class
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS%in%"TRAIN")]
	qFile <- sprintf("%s/%s_query",pDir,g)
	GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
	
	resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)

	predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
}

## ----label-predictions, eval=TRUE,cache=FALSE----------------------------
predClass <- GM_OneVAll_getClass(predRes)
cat("Predicted classes\n")

## ----eval-perf, eval=TRUE,cache=FALSE------------------------------------
both <- merge(x=pheno,y=predClass,by="ID")
print(table(both[,c("STATUS","PRED_CLASS")]))
pos <- (both$STATUS %in% "LumA")
tp <- sum(both$PRED_CLASS[pos]=="LumA")
fp <- sum(both$PRED_CLASS[!pos]=="LumA")
tn <- sum(both$PRED_CLASS[!pos]=="other")
fn <- sum(both$PRED_CLASS[pos]=="other")
cat(sprintf("Accuracy = %i of %i (%i %%)\n",tp+tn,nrow(both),
			round(((tp+tn)/nrow(both))*100)))
cat(sprintf("PPV = %i %%\n", round((tp/(tp+fp))*100)))
cat(sprintf("Recall = %i %%\n", round((tp/(tp+fn))*100)))

}, error=function(ex) {
	print(ex)
}, finally={
	cat("Closing log.\n")
	sink(NULL)
})
