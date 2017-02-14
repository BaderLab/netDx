#' Build predictor with n-way resampling, for CNV-based nets
#'
#' @param pheno (data.frame) patient IDs (ID) and class of interest (STATUS)
#' @param p_GR (GRanges) patient CNVs. Rows are patients, columns are unit
#' measures (e.g. genes)
#' @param predClass (char) class to build predictor for. Must be a value
#' in the STATUS column of the pheno param
#' @param unitSet_GR (list) unit groupings, each of which will be 
#' converted to its own patient similarity network. 
#' @param unit_GR (GRanges) genomic locations of units (e.g. genes).
#' Provided data in this form speeds test evaluation
#' @param pctT (numeric 0 to 1) fraction of samples for training
#' @param numResamples (integer) number of resamplings for cross validation
#' @param cliqueReps (integer) number of permutations to perform in
#' clique filtering
#' @param nFoldCV (integer) number of folds for cross-validation within each
#' resampling
#' @param filter_WtSum (numeric between 5-100) Limit to top-ranked 
#' networks such that cumulative weight is less than this parameter. 
#' e.g. If filter_WtSum=20, first order networks by decreasing weight; 
#' then keep those whose cumulative weight <= 20.
#' @param numCores (integer) number of parallel processing jobs
#' @param GMmemory (integer) memory (Gb) for GeneMANIA to run.
#' @param outDir (char) path to output directory
#' @param overwrite (logical) Overwrite existing results in outDir? By 
#' default, no, return warning; but set to TRUE if outDir should be deleted
#' and populated with fresh results from this function.
#' @param seed_trainTest (integer) RNG seed for initial train/test split.
#' Set to NULL to use current state of RNG. Set to a constant for 
#' reproducibility
#' @param seed_resampling (integer) RNG seed for training assignment for
#' of the resamplings (splitTrainTest_partition() method)
#' @param seed_CVqueries (integer) RNG seed for queries in k-fold cross
#' validation (used by makeCVqueries() which is called by 
#' GM_runCV_featureSet())
#' @param ... parameters for Nway_netSum()
#' @return (list) Includes:
#' * pheno: (data.frame) phenotype matrix as used by predictor
#' * netScores: (data.frame) networks and scores following feature 
#' selection. Nets with score>0 are shown.
#' * resamplingResults (list) performance of resampling rounds. Output of
#' resampling_pickBestCutoff_CNV()
#' * cliqueNets (char) vector of nets that pass clique filtering in 1+
#' round of resampling
#' * testRes (list) results of evaluating performance on the test partition.
#' output of resampling_predTest_CNV()
#' @export
buildPredictor_resampling_CNV <- function(pheno, p_GR, predClass, 
    unitSet_GR,unit_GR,pctT=0.7,numResamples=3L, cliqueReps=500L,
	nFoldCV=10L,filter_WtSum=100L, numCores=1L,GMmemory=4L,
	outDir=".",overwrite=FALSE,seed_trainTest=42L,seed_resampling=42L,
	seed_CVqueries=42L,...) {
	
TEST_MODE <- FALSE

if (!TEST_MODE) {
if (file.exists(outDir)) {
	if (!overwrite) stop("output directory exists. Choices: provide path to non-existing directory, set overwrite option to TRUE, or manually delete this directory")
	else unlink(outDir,recursive=TRUE)
}
dir.create(outDir)
}

pheno$STATUS[which(!pheno$STATUS %in% predClass)] <- "other"
subtypes <- c(predClass,"other")

if (pctT < 1) {
pheno$TT_STATUS <- splitTestTrain(pheno,
    	pctT=pctT,setSeed=seed_trainTest)
} else {
	warning("No test partition. Entire dataset will be used to compute test scores");
	pheno$TT_STATUS <- "TRAIN"
}

pheno_FULL	<- pheno
pGR_FULL 	<- p_GR
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
p_GR		<- p_GR[which(p_GR$ID %in% pheno$ID)]

cat("Training samples\n")
print(table(pheno$STATUS))

if (!TEST_MODE){
# --------------------------------------------------
# Phase 1. Feature selection, assigning net scores
# --------------------------------------------------
cat("* Make CNV-based patient networks\n")
netDir <- sprintf("%s/networks_orig",outDir)
print(numCores)
netList <- makePSN_RangeSets(p_GR, unitSet_GR,netDir,verbose=FALSE,
	 numCores=numCores)

cat("* Counting patients in net\n")
t0 <- Sys.time()
p 	<- countPatientsInNet(netDir,netList, pheno$ID)
t1 <- Sys.time()
print(t1-t0)
tmp	<- updateNets(p,pheno,writeNewNets=FALSE)
p		<- tmp[[1]]
pheno	<- tmp[[2]] 

pheno_train_in_nets <- pheno
dat <- list(
	pheno_FULL=pheno_FULL,
	pheno_training_in_nets=pheno,
	netmat_training=p)
save(dat,file=sprintf("%s/pheno_main.Rdata",outDir))
rm(dat)

# compute net scores over resampled training data
cat("* Running N-way resampling for feature scores\n")
t0 <- Sys.time()
Nway_netSum(p,pheno,predClass=predClass,outDir,netDir,
			cliqueReps=cliqueReps,numCores=numCores,
			splitN=numResamples,seed_resampling=seed_resampling,
			nFoldCV=nFoldCV,filter_WtSum=filter_WtSum,
			seed_CVqueries=seed_CVqueries,
			...)
t1 <- Sys.time()
print(t1-t0)
} else {
	warning("TEST_MODE enabled: not computing net scores\n")
}

# --------------------------------------------------
# Phase 2. Get cutoff with best accuracy
# --------------------------------------------------

# note that the patients overlapping called nets has already been calculated by
# nWay_netSum::RR_featureTally(). We are reading this output and computing
# accuracy for each cutoff
#
# TODO - nWay_netSum should just pass the resampPerf data as output
# instead of this function having to load an output file.
# pick best cutoff

cat("* Picking best cutoff\n")
cat("------------------------------------\n")
load(sprintf("%s/resamplingPerf.Rdata",outDir))
pdf(sprintf("%s/train_accuracy.pdf",outDir))
tryCatch({
    res <- resampling_pickBestCutoff_CNV(resampPerf[["cliqueNets"]])
},error=function(ex) {
    print(ex)
},finally={
    dev.off()
})
resamplingRes <- res
cat(sprintf("Best cutoff is: %i\n", res$bestCutoff))

# load feature-selected net names
# TODO this should be provided by nWay_netSum
pTally <- read.delim(sprintf("%s/pathway_cumTally.txt",outDir),
	sep="\t",h=T,as.is=T)
pTally_full <- pTally

# fetch cliqueNets
# TODO this should be provided by nWay_netSum
e1 <- new.env() # must load in a contained environment to avoid
				# potential variable overwrite
load(sprintf("%s/resampling_savedDat.Rdata",outDir))
cliqueNets <- get('cliqueNets', e1)
rm(e1)
cliqueNets <- unique(unlist(cliqueNets))
cliqueNets <- sub("_cont.txt","",cliqueNets)

if (pctT == 1) {
	cat("No TEST partition. Returning results\n")
	out <- list(
		pheno=pheno_FULL,
		netScores=pTally_full,
		resamplingResults=resamplingRes,
	##	bestCutoff=bestCutoff,
		cliqueNets=cliqueNets
	)
	return(out)
}

# --------------------------------------------------
# Phase 3. Compute accuracy for test patients
# --------------------------------------------------
pheno <- subset(pheno_FULL, TT_STATUS %in% "TEST")
p_GR <- pGR_FULL[which(pGR_FULL$ID %in% pheno$ID)]

cat(sprintf("Test: %i patients ; %i ranges\n", nrow(pheno),length(p_GR)))
cat("Class breakdown:\n")
print(table(pheno[,c("STATUS")]))

testRes <- list() # performance of test for each cutoff

for (curr_cutoff in 1:max(pTally[,2])) {
	cat(sprintf("Test: score=%i\n", curr_cutoff))
	pTally <- pTally_full
	pTally <- pTally[which(pTally[,2]>=curr_cutoff),1]
	pTally <- sub("_cont.txt","",pTally)
	cat(sprintf("\t%i pathways\n", length(pTally)))

	idx <- which(names(unitSet_GR)%in% cliqueNets)
	if (length(idx)<length(cliqueNets)) {
    	cat("some cliqueNets not located\n")
    	browser()
	}
	cNetGenes <- lapply(unitSet_GR[idx],function(x) { x$name })
	cNetGenes <- unit_GR[which(unit_GR$name %in% unlist(cNetGenes))]
	cat(sprintf("\t%i units in clique-filtered nets\n",length(cNetGenes)))

	# fetch selFeature_GR
	# now also get genes in feature-selected pathways
	idx <- which(names(unitSet_GR) %in% pTally)
	if (length(idx) < length(pTally)) {
   	 cat("some pTally not located\n")
   	 browser()
	}	
	selGR_genes <- lapply(unitSet_GR[idx], function(x) {x$name})
	selGR_genes <- unit_GR[which(unit_GR$name %in% unlist(selGR_genes))]
	cat(sprintf("\t%i units in feature-selected nets\n",length(selGR_genes))) 

	# count predClass/other overlapping feature selected pathways
	testRes[[curr_cutoff]] <- resampling_predTest_CNV(
		pheno,p_GR,
		selFeature_GR=selGR_genes,predClass=predClass,
    	cliqueNetGenes=cNetGenes)
##	bestCutoff <- resamplingRes$bestCutoff

	#cat(sprintf("Test performance: Cutoff=%i; Acc=%1.2f ; PPV=%1.2f\n",
	#	bestCutoff, testRes$acc, testRes$ppv))
	rm(cNetGenes,selGR_genes,pTally)
}
	save(testRes,file="testPerformance.Rdata")
	save(testRes,pheno,p_GR,cliqueNets,
	file=sprintf("%s/testPerformance.Rdata",outDir))

out <- list(
	pheno=pheno_FULL,
	netScores=pTally_full,
	resamplingResults=resamplingRes,
##	bestCutoff=bestCutoff,
	cliqueNets=cliqueNets,
	testRes=testRes
)

out
}
