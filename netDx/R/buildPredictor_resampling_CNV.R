#' Build predictor with n-way resampling, for CNV-based nets
#'
#' @param pheno (data.frame) patient IDs (ID) and class of interest (STATUS)
#' @param p_GR (GRanges) patient CNVs. Rows are patients, columns are unit
#' measures (e.g. genes)
#' @param predClass (char) class to build predictor for. Must be a value
#' in the STATUS column of the pheno param
#' @param unitSet_GR (list) unit groupings, each of which will be 
#' converted to its own patient similarity network. 
#' @param pctT (numeric 0 to 1) fraction of samples for training
#' @param numResamples (integer) number of resamplings for cross validation
#' @param cliqueReps (integer) number of permutations to perform in
#' clique filtering
#' @param nFoldCV (integer) number of folds for cross-validation within each
#' resampling
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
#' * TT_STATUS: initial train/test split assignment
#' * netScores: (list) keys are classes, value is a data.frame with 
#' networks and scores following feature selection. Nets with score>0 are
#' shown.
#' * confmat: (matrix) confusion matrix with one row per cutoff
#' * predRanks: (list) GM rankings of test samples for each network cutoff
#' In addition, all the intervening work will be stored in <outDir>
#' @export
buildPredictor_resampling_CNV <- function(pheno, p_GR, predClass, 
    unitSet_GR,pctT=0.7,numResamples=3L, cliqueReps=500L,
	nFoldCV=10L,numCores=1L,GMmemory=4L,
	outDir=".",overwrite=FALSE,seed_trainTest=42L,seed_resampling=103L,
	seed_CVqueries=42L,...) {
	
if (file.exists(outDir)) {
	if (!overwrite) stop("output directory exists. Choices: provide path to non-existing directory, set overwrite option to TRUE, or manually delete this directory")
	else unlink(outDir,recursive=TRUE)
}
dir.create(outDir)

pheno$STATUS[which(!pheno$STATUS %in% predClass)] <- "other"
subtypes <- c(predClass,"other")

pheno$TT_STATUS <- splitTestTrain(pheno,
    	pctT=pctT,setSeed=seed_trainTest,predClass=predClass)

pheno_FULL	<- pheno
pGR_FULL 	<- p_GR
pheno		<- subset(pheno,TT_STATUS %in% "TRAIN")
pGR_FULL	<- p_GR[which(p_GR$ID %in% pheno$ID)]

cat("Training samples\n")
print(table(pheno$STATUS))

# --------------------------------------------------
# Phase 1. Feature selection, assigning net scores
# --------------------------------------------------
cat("* Make CNV-based patient networks\n")
netDir <- sprintf("%s/networks_orig",outDir)
print(numCores)
netList <- makePSN_RangeSets(p_GR, unitSet_GR,netDir,verbose=FALSE,
							 numCores=numCores)

p 	<- countPatientsInNet(netDir,netList, pheno$ID)
tmp	<- updateNets(p,pheno,writeNewNets=FALSE)
p		<- tmp[[1]]
pheno	<- tmp[[2]] 

Nway_netSum(p,pheno,predClass=predictClass,outDir,netDir,
			cliqueReps=cliqueReps,numCores=numCores,
			splitN=numResamples,nFoldCV=nFoldCV,...)

browser()

}
