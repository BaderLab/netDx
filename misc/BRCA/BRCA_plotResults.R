#' plot BRCA results

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_BRCA/output/BRCA_180117"

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

xpr_genes <- rownames(xpr)
pathwayList <- lapply(pathwayList,function(x) x[which(x %in% xpr_genes)])

pheno$STATUS[which(pheno$STATUS!="LumA")] <- "other"

#out <- plotAllResults(pheno, sprintf("%s/pred",rootDir),
#							 outDir=sprintf("%s/plot",rootDir),
#               fsCutoff=10,fsPctPass=0.7,pathwaySet=pathwayList)

inDir <- sprintf("%s/pred",rootDir)
outDir <- sprintf("%s/plot",rootDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()
featScores <- getFeatureScores(inDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
})

netInfoFile <- sprintf("%s/inputNets.txt",inDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir)
