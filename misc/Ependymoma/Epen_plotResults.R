#' plot BRCA results
require(netDx)
require(netDx.examples)
rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_Ependymoma"
dataDir <- sprintf("%s/input/data_for_shraddha/original_data/Toronto-comparison2-without-spinals",rootDir)

xpr <- read.delim(sprintf("%s/TOR-ST-PFPURE-PFMIX-SEP16.gct",dataDir),skip=2,h=T,as.is=T)
rownames(xpr) <- xpr[,1]
xpr <- xpr[,-(1:2)]
sampType <- scan(sprintf("%s/TOR-ST-PFPURE-PFMIX-SEP16.cls",dataDir),skip=2)
sampType <- as.integer(sampType)
# from Ruth
#  st = 0, PFPURE = 1 and PFMIX = 2
pheno <- data.frame(ID=colnames(xpr),INT_STATUS=sampType)
pheno$ID <- as.character(pheno$ID)
st <- c("ST","PFPURE","PFMIX")
pheno$STATUS <- st[sampType+1]

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

inDir <- sprintf("%s/output/Epen_180118/pred",rootDir)
outDir <- sprintf("%s/output/Epen_180118/plot",rootDir)
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
