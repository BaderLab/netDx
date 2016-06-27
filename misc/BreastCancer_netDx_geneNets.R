## netDx BreastCancer. gene expression only.
## gene-level networks.

## ----set-params,cache=FALSE,eval=TRUE------------------------------------
rm(list=ls())

# Change this to a local directory where you have write permission
outDir <- "~/tmp/TCGA_BRCA_geneXpr_geneNets" 

numCores 	<- 8L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
cutoff		<- 9L  	# score cutoff for feature-selected networks
TRAIN_PROP <- 0.67 	# fraction of samples to use for training

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

## ----load-packages, cache=FALSE,eval=TRUE--------------------------------
require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

sink(sprintf("%s/BreastCancer_GeneExprOnly.log",outDir),split=TRUE)
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

nr <- nrow(xpr)
cat("*** Excluding missing data\n")
xpr <- na.omit(xpr)
cat(sprintf("%i of %i genes excluded\n",nr-nrow(xpr), nr))

# fix gene name with slash
idx <- grep("NOP5/NOP58",rownames(xpr))
rownames(xpr)[idx] <- "NOP58"

## ----read-pathways, cache=FALSE,eval=TRUE--------------------------------
cat("* Creating gene-sets, each with one gene\n")
pathwayList <- as.list(rownames(xpr))
names(pathwayList) <- rownames(xpr)

head(pathwayList)

## ----make-xpr-psn, eval=TRUE---------------------------------------------
profDir <- sprintf("%s/profiles",outDir)
netDir <- sprintf("%s/networks",outDir)

geneSim <- function(x) {
    if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x)-min(x)
    
    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(g[i]-g[j])/(max(g)-min(g)))
    # where g is the eMB.xpression vector for each gene
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}

netList <- makePSN_NamedMatrix(xpr, rownames(xpr), 
        pathwayList,profDir,verbose=TRUE,
        numCores=numCores,
		simMetric="custom", customFunc=geneSim,
		sparsify=TRUE)
netList <- unlist(netList)
head(netList)

## ----make-cnv-psn, cache=FALSE,eval=TRUE---------------------------------
data(genes)
gene_GR     <- GRanges(genes$chrom,
   IRanges(genes$txStart,genes$txEnd),
   name=genes$name2)
###path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)

###names(path_GRList) <- paste("CNV_",names(path_GRList),sep="")
## warning: this step can take 2-5 minutes depending on the
## number of processes running in parallel
###netList2 <- makePSN_RangeSets(cnv_GR, path_GRList,profDir,verbose=F)
###cat(sprintf("CNV: Got %i networks\n",length(netList2)))
###
##### ----cnv-psn-look, cache=FALSE,eval=TRUE---------------------------------
###head(unlist(netList2))

## ----create-gm-db, eval=TRUE---------------------------------------------
# now create database
dbDir	<- GM_createDB(profDir, pheno$ID, outDir,numCores=numCores)

## ----feature-selection,cache=FALSE,eval=TRUE-----------------------------
## repeat process for each class
for (g in subtypes) {
    pDir <- sprintf("%s/%s",outDir,g)
    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
	dir.create(pDir)

	cat(sprintf("\n******\nSubtype %s\n",g))
	pheno_subtype <- pheno
	
	## label patients not in the current class as a residual
	pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "nonpred"
	## sanity check
	print(table(pheno_subtype$STATUS,useNA="always"))
    
	resDir    <- sprintf("%s/GM_results",pDir)
	## query for feature selection comprises of training 
	## samples from the class of interest
	trainPred <- pheno$ID[which(pheno$STATUS %in% g)]
	
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
		profDir,verbose=F,numCores=numCores,
		simMetric="custom",customFunc=geneSim, 
		sparsify=TRUE)
	###tmp <- makePSN_RangeSets(cnv_FULL,
	###	path_GRList[which(names(path_GRList)%in% pTally)],
	###		profDir,verbose=FALSE)
	# create db
	dbDir <- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

	# query of all training samples for this class
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
							 pheno$TT_STATUS%in%"TRAIN")]
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
