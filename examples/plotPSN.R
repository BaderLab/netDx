suppressMessages(require(netDx))
suppressMessages(require(netDx.examples))

phenoFile <- sprintf("%s/extdata/KIRC_pheno.rda",
	path.package("netDx.examples"))
lnames <- load(phenoFile)
outDir <- paste(getwd(),"plots",sep="/")
if (!file.exists(outDir)) dir.create(outDir)
#setwd(outDir)

inDir <- sprintf("%s/extdata/KIRC_output",
	path.package("netDx.examples"))
all_rngs <- list.dirs(inDir, recursive = FALSE)

predClasses <- c("SURVIVEYES","SURVIVENO")
predFiles <- unlist(lapply(all_rngs, function(x) 
		paste(x, "predictionResults.txt", sep = "/")))
featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))

featSelNet <- lapply(featScores, function(x) {
	callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
})

netInfo <- plotIntegratedPSN(pheno=pheno,baseDir=sprintf("%s/rng1",inDir),
	netNames=featSelNet,outDir=outDir,edgeStroke="#000055",
	edgeTransparency=140,nodePal="Spectral",nodeSize=200,
	nodeTransparency=200,edgeWidth=2,imageFormat="jpeg")
print(netInfo)

sessionInfo()

