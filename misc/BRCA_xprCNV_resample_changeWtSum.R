#' LumA binary classification with resampling.
#' Work which loops over different levels of network exclusion in inner
#' cross validation loop.
#' 
#' Set to run multiple times on scinet via command-line args
#' Arguments:
#' 1. (char) outDir path to output directory
#' 2. (integer) seed for train/test split
#' 3. (integer) seed for resampling split 
args    <- commandArgs(TRUE)
print(args)

DEBUG_MODE <- FALSE

outDir          <- args[1]
seed_trainTest  <- as.integer(args[2])
seed_resampling <- as.integer(args[3])
filter_WtSum	<- as.integer(args[4])
#outDir          <- args[1]
#seed_trainTest  <- seed1
#seed_resampling <- seed2

# ---------------------------------------------------
# do work

#if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
#dir.create(outDir)

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
	cat(sprintf("RNG seeds: Train-test split = %i ; Resampling = %i\n",
		seed_trainTest, seed_resampling))
	cat(sprintf("Filter wt sum = %i\n", filter_WtSum))

	pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
 	   path.package("netDx.examples"))
	pathwayList <- readPathways(pathFile)

	data(genes)
	gene_GR     <- GRanges(genes$chrom,
  						IRanges(genes$txStart,genes$txEnd),
   						name=genes$name2)

	cat("* Limiting to pathway genes\n")
	path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
	names(path_GRList) <- paste("CNV_",names(path_GRList),sep="")

	predDir <- sprintf("%s/predictor",outDir)

	if (DEBUG_MODE) {
		warning("** in debug mode **")
		Sys.sleep(2)
		nFold <- 2L
		numResamp <- 2L
		num <- round(0.2*length(pathwayList))
		pathwayList <- pathwayList[1:num]
	} else {
		nFold <- 10L
		numResamp <- 3L
	}

	cat("* Mapping CNV to genes\n")
	cnv_GR <- getRegionOL(cnv_GR, path_GRList)

	cat("* Running predictor\n")
	t0 <- Sys.time()
	out <- buildPredictor_resampling_Mixed(pheno=ph,pdat=xpr,
		predClass="LumA",
		p_GR=cnv_GR, unitSet_GR=path_GRList,
		nFoldCV=nFold, numResamples=numResamp,
		unitSets=pathwayList,numCores=numCores,outDir=predDir,
		overwrite=TRUE,GMmemory=GMmemory,
		seed_trainTest=seed_trainTest, seed_resampling=seed_resampling,
		filter_WtSum=filter_WtSum)
	save(out,file=sprintf("%s/FinalResults.Rdata",outDir))
	print(Sys.time()-t0)

}, error=function(ex) {
	print(ex)
}, finally={
	cat("Closing log.\n")
	print(Sys.time())
	sink(NULL)
})
