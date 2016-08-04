#' LumA binary classification with resampling.
#' Set to run multiple times on scinet via command-line args
#' Arguments:
#' 1. (char) outDir path to output directory
#' 2. (integer) seed for train/test split
#' 3. (integer) seed for resampling split 
args    <- commandArgs(TRUE)
# Change this to a local directory where you have write permission

#runMe <- function(outDir,seed1,seed2) {
outDir          <- args[1]
seed_trainTest  <- as.integer(args[2])
seed_resampling <- as.integer(args[3])
#outDir          <- args[1]
#seed_trainTest  <- seed1
#seed_resampling <- seed2

# ---------------------------------------------------
# do work

if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

numCores 	<- 8L  	# num cores available for parallel processing
GMmemory 	<- 4L  	# java memory in Gb
trainProp	<- 0.67 # fraction of samples to use for training

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)
ph <- pheno; rm(pheno)

sink(sprintf("%s/BreastCancer_GeneExprOnly.log",outDir),split=TRUE)
tryCatch({
	cat("Start time\n")
	print(Sys.time())

	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
 	   path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)

	### TODO Should be able to pass GM memory setting.
	out <- buildPredictor_resampling(pheno=ph,pdat=xpr,predClass="LumA",
		nFoldCV=10L, numResamples=3L,
		unitSets=pathwayList,numCores=8L,outDir=outDir,overwrite=TRUE,
		seed_trainTest=seed_trainTest, seed_resampling=seed_resampling)
	save(out,file=sprintf("%s/FinalResults.Rdata",outDir))

    ### now clean up intermediate files 
    system(sprintf("rm -r %s/LumA %s/other %s/eval %s/test %s/dataset %s/tmp",
		outDir,outDir,outDir,outDir,outDir,outDir))
    

}, error=function(ex) {
	print(ex)
}, finally={
	cat("Closing log.\n")
	print(Sys.time())
	sink(NULL)
})
}
