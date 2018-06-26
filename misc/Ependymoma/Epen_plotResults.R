#' plot BRCA results
rm(list=ls())
require(netDx)
require(netDx.examples)
#rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_Ependymoma"
rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"

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

#setName <- "Epen_prunedOneNet_0.001_180409"
#setName <- "Epen_prunedPathway_0.01_180409"
#setName <- "Epen_prunedOneNet_180522"
#setName <- "Epen_prunedOneNet_0.001_180523"
#setName <- "Epen_nopruneOneNet_180523"
setName <- "Epen_lassoOutsideCV_180523"
inDir <- sprintf("%s/output/%s/pred",rootDir,setName)
outDir <- sprintf("%s/output/%s/plot",rootDir,setName)

if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()
auroc <- unlist(lapply(predPerf, function(x) x$auroc))
aupr <- unlist(lapply(predPerf, function(x) x$aupr))
acc <- unlist(lapply(predPerf, function(x) x$accuracy))

den <- sqrt(length(auroc))
cat("--------------\n")
cat(sprintf("Performance: %s\n",setName))
cat(sprintf("AUROC = %1.2f +/- %1.2f\n",mean(auroc),sd(auroc)/den))
cat(sprintf("AUPR = %1.2f +/- %1.2f\n",mean(aupr),sd(aupr)/den))
cat(sprintf("Accuracy = %1.2f +/- %1.2f\n",mean(acc),sd(acc)/den))
cat("--------------\n")


###featScores <- getFeatureScores(inDir,predClasses=predClasses)
###featSelNet <- lapply(featScores, function(x) {
###    callFeatSel(x, fsCutoff=10, fsPctPass=0.9)
###})
###
###netInfoFile <- sprintf("%s/inputNets.txt",inDir)
###netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
###EMap_input <- writeEMapInput_many(featScores,pathwayList,
###      netInfo,outDir=outDir)
