#' scratchpad for fixing predictor builder.
rm(list=ls())

outDir <- "/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test"
numResamples <- 2L #L
nFoldCV <- 10L
subtypes <- c("LumA","other")
numCores <- 8L
predClass <- "LumA"

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
   path.package("netDx.examples"))
unitSets <- readPathways(pathFile)

maxScore <- 3 ###numResamples * nFoldCV

## load TT_status etc
load(sprintf("%s/predResults.Rdata",outDir))
pheno_FULL <- out$pheno
xpr_FULL	<- xpr
pheno_train <- subset(pheno_FULL, TT_STATUS %in% "TRAIN")
xpr_train <- xpr[,which(colnames(xpr)%in% pheno_train$ID)]

##  pathways by subtype
TT_resamp <- list() # train/test status for resampling
netScores <- list()
for (g in subtypes) {
	pTally <- read.delim(
		sprintf("%s/%s/%s_pathwayScore.txt",outDir,g,g),
		sep="\t",h=T,as.is=T)
	pTally <- pTally[which(pTally[,2]>=1),]
	pTally[,1] <- sub(".profile","",pTally[,1])
	netScores[[g]] <- pTally

	load(sprintf("%s/%s/TT_STATUS_resampling.Rdata",outDir,g))
	TT_resamp[[g]] <- TT_STATUS
}

###perf_resamp <- list()
###for (rep in 1:numResamples) {
###	## set pheno so that query is limited to samples that were training
###	## for resampling i
###	## make sure that pdat is limited to samples in pheno (Should be all
###	## training samples)
###	perfDir <- sprintf("%s/eval/part%i",outDir,rep)
###	if (file.exists(perfDir)) unlink(perfDir,recursive=TRUE)
###	dir.create(perfDir,recursive=TRUE)
###	pheno_cur <- pheno_train
###	pheno_cur$TT_STATUS <- TT_resamp[[g]][[rep]]
###	print(table(pheno_cur[,c("STATUS","TT_STATUS")]))
###	perf_resamp[[rep]] <- GM_predClass_cutoffs(pheno_cur,xpr_train, 
###		predClass=predClass,netScores=netScores,unitSets=unitSets,
###		maxScore=maxScore,outDir=perfDir,numCores=numCores)
###}
###save(perf_resamp,file=sprintf("%s/resamplingPerf.Rdata",outDir))
load(sprintf("%s/resamplingPerf.Rdata",outDir))

# pick the cutoff that has the best mean performance across resamplings
perfmat <- matrix(NA, nrow=nrow(perf_resamp[[1]][["confmat"]]),
	ncol=4*length(perf_resamp))
colnames(perfmat) <- rep(c("tpr","fpr","accuracy","ppv"),
	 length(perf_resamp))
ctr <- 1
for (k in 1:length(perf_resamp)) {
	perfmat[,ctr:(ctr+3)] <- perf_resamp[[k]][["confmat"]][,6:9]
	ctr <- ctr+4
}
mu <- matrix(NA,nrow=nrow(perfmat),ncol=4)
for (k in 1:4) {
	print(seq(k,ncol(perfmat),4))
	mu[,k] <- rowMeans(perfmat[,seq(k,ncol(perfmat),4)])
}
colnames(mu) <- c("tpr","fpr","accuracy","ppv")

bestIdx <- which.max(mu[,3]) # cutoff is the one with best accuracy
bestCutoff <- perf_resamp[[1]][["confmat"]][bestIdx,1]
cat(sprintf("Best cutoff = %i ; TPR=%1.2f; FPR=%1.2f; Accuracy=%1.2f; PPV=%1.2f\n",
	bestCutoff, mu[bestIdx,1],mu[bestIdx,2],mu[bestIdx,3],mu[bestIdx,4]))	

cat("\n\n* Model evaluation\n")
finalDir <- sprintf("%s/test", outDir)
finalNets <- lapply(netScores, function(x) x[which(x[,2]>=bestCutoff),1])
testRes <- GM_predClass_once(pheno_FULL,xpr_FULL,predClass,unitSets,
	patNets=finalNets,outDir=finalDir,numCores=numCores)





