#' TCGA Breast Cancer use-case two-class classification
#' 
#' @details Class of interest: Luminal A. Data types are gene expression 
#' and CNV
rm(list=ls())

# Change this to a local directory where you have write permission
outDir <- "~/tmp/TCGA_BRCA" 

numCores 	<- 1L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
cutoff		<- 9L  	# score cutoff for feature-selected networks
TRAIN_PROP <- 0.67 	# fraction of samples to use for training

require(netDx)
require(netDx.examples)

data(TCGA_BRCA)

# --------------------------------------------------------
# begin patient networks work
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

subtypes<- c("LumA")

pheno$STATUS[which(!pheno$STATUS %in% subtypes)] <- "other"
mega_TT <- rep("TRAIN",nrow(pheno))
# first separate the test samples
for (uq in unique(pheno$STATUS)) {
	idx <- which(pheno$STATUS %in% uq)
	vec <- rep("TRAIN",length(idx))
	k <- floor((1-TRAIN_PROP)*length(idx))
	vec[sample(1:length(idx),k,F)]<- "TEST"
	mega_TT[idx] <- vec
}
cat("Train/test splits\n")
pheno$TT_STATUS <- mega_TT
print(table(pheno[,c("STATUS","TT_STATUS")]))
subtypes <- c(subtypes,"other") # add residual

pheno_FULL	<- pheno
xpr_FULL 	<- xpr
cnv_FULL	<- cnv_GR
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
xpr			<- xpr[,which(colnames(xpr)%in% pheno$ID)]
cnv_GR		<- cnv_GR[which(cnv_GR$ID %in% pheno$ID)]

cat(sprintf("Limiting to TRAIN patients: %i of %i\n",nrow(pheno),
			nrow(pheno_FULL)))

cat("\tMake patient networks\n")
# make pathway list
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
	path.package("netDx.examples"))
pathwayList    <- readPathways(pathFile)

# create patient networks
profDir <- sprintf("%s/profiles",outDir)
netDir <- sprintf("%s/networks",outDir)
## compare num interactions of a network by own corr vs GeneMANIA corr.

netList <- unlist(makePSN_NamedMatrix(xpr, rownames(xpr), pathwayList,
		profDir,verbose=FALSE,numCores=numCores,writeProfiles=TRUE))

# add CNV nets
data(genes)
gene_GR     <- GRanges(genes$chrom,
					   IRanges(genes$txStart,genes$txEnd),
                       name=genes$name2)
path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
names(path_GRList) <- paste("CNV_",names(path_GRList),sep="")
netList2 <- makePSN_RangeSets(cnv_GR, path_GRList,profDir,verbose=F)
cat(sprintf("CNV: Got %i networks\n",length(netList2)))

# now create database
dbDir	<- GM_createDB(profDir, pheno$ID, outDir,numCores=numCores)

for (g in subtypes) {
    pDir <- sprintf("%s/%s",outDir,g)
    if (file.exists(pDir)) unlink(pDir,recursive=TRUE)
	dir.create(pDir)

	cat(sprintf("Subtype %s\n",g))
	pheno_subtype <- pheno
	pheno_subtype$STATUS[which(!pheno_subtype$STATUS %in% g)] <- "other"
	print(table(pheno_subtype$STATUS,useNA="always"))
    
	resDir    <- sprintf("%s/GM_results",pDir)
	# samples from which we want to make queries for cross-validation
	# training samples for predictor class
	trainPred <- pheno$ID[which(pheno$STATUS %in% g)]
	GM_runCV_featureSet(trainPred, resDir, dbDir$dbDir, 
		nrow(pheno_subtype),verbose=T, numCores=numCores,
		GMmemory=GMmemory)
	
	prank <- dir(path=resDir,pattern="PRANK$")
	nrank <- dir(path=resDir,pattern="NRANK$")
	cat(sprintf("Got %i prank files\n",length(prank)))
	GMres <- list()
	for (f in prank)  {
		GMres[[f]] <- GM_getQueryROC(sprintf("%s/%s",resDir,f),
									   pheno_subtype,g)
	}
	outFile <- sprintf("%s/GM_perf.Rdata",pDir)
	save(GMres,file=outFile)
		rm(GMres)
	pTally		<- GM_networkTally(paste(resDir,nrank,sep="/"))
	tallyFile	<- sprintf("%s/%s_pathway_CV_score.txt",resDir,g)
	write.table(pTally,file=tallyFile,sep="\t",col=T,row=F,quote=F)
	fs<- pTally[which(pTally[,2]>=cutoff),]
	fs[,1]<-sub(".profile","",fs[,1])
	fs[,1] <- sub("_cont","",fs[,1])

	outFile <-sprintf("%s/%s_cutoff%i_nets.gmt",resDir,g,cutoff)
	system(sprintf("cat /dev/null > %s",outFile))
	for (k in fs[,1]) {
		if (any(grep("^CNV_",k))) {
			m <- sub("CNV_","",k)
		} else {
			m <- k
		}
		#only keep genes present in our data
		g <- intersect(pathwayList[[m]],rownames(xpr))
		cat(sprintf("%s\t%s\t%s\n",k,k,paste(g,collapse="\t")),
			file=outFile,append=TRUE)
	}
}

# now create GM databases for each type
# should contain train + test patients
pheno <- pheno_FULL
predRes <- list()
for (g in subtypes) {
	pDir <- sprintf("%s/%s",outDir,g)
	pTally <- read.delim(
		sprintf("%s/GM_results/%s_pathway_CV_score.txt",pDir,g),
		sep="\t",h=T,as.is=T)
	pTally <- pTally[which(pTally[,2]>=cutoff),1]
	pTally <- sub(".profile","",pTally)
	pTally <- sub("_cont","",pTally)

	cat(sprintf("%s: %i pathways\n",g,length(pTally)))
	profDir <- sprintf("%s/profiles",pDir)

	# make nets
	netList <- unlist(makePSN_NamedMatrix(xpr_FULL,rownames(xpr),
		pathwayList[which(names(pathwayList)%in% pTally)],
		profDir,verbose=F,numCores=numCores,writeProfiles=TRUE))
	netList2<- makePSN_RangeSets(cnv_FULL,
		path_GRList[which(names(path_GRList)%in% pTally)],
		profDir,verbose=FALSE)
	# create db
	dbDir <- GM_createDB(profDir,pheno$ID,pDir,numCores=numCores)

	# query is all training samples
	qSamps <- pheno$ID[which(pheno$STATUS %in% g & 
							 pheno$TT_STATUS%in%"TRAIN")]
	cat(sprintf("%s: %i items in query\n",g,length(qSamps)))	 
	qFile <- sprintf("%s/%s_query",pDir,g)
	GM_writeQueryFile(qSamps,"all",nrow(pheno),qFile)
	
	resFile <- runGeneMANIA(dbDir$dbDir,qFile,resDir=pDir)

	predRes[[g]] <- GM_getQueryROC(sprintf("%s.PRANK",resFile),pheno,g)
	
}
save(predRes,file=sprintf("%s/testPredictions.Rdata",outDir))

cat("------------------\nPrediction stats\n")
# now rank each patient
predClass <- GM_OneVAll_getClass(predRes)
cat("Predicted classes\n")
print(table(predClass[,c("PRED_CLASS")]))
write.table(predClass,file=sprintf("%s/predClass_all.txt",outDir),sep="\t",
			col=TRUE,row=F,quote=F)

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

