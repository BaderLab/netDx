#' plot BRCA results
rm(list=ls())
require(netDx)
require(netDx.examples)
rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_Ependymoma"

inFile <- sprintf("%s/input/netDx_prepared/Ependymoma_cohortMerged_180125.Rdata",rootDir)
load(inFile)

# exclude ST
idx <- which(pheno$STATUS=="ST") 
pheno <- pheno[-idx,]
xpr <- xpr[,-idx]

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

#out <- plotAllResults(pheno, sprintf("%s/pred",rootDir),
#							 outDir=sprintf("%s/plot",rootDir),
#               fsCutoff=10,fsPctPass=0.7,pathwaySet=pathwayList)

inDir <- sprintf("%s/output/Epen_180125/pred",rootDir)
outDir <- sprintf("%s/output/Epen_180125/plot",rootDir)
if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()
featScores <- getFeatureScores(inDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.9)
})

netInfoFile <- sprintf("%s/inputNets.txt",inDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir)
