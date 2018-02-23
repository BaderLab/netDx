#' plot BRCA results
rm(list=ls())
require(netDx)
require(netDx.examples)
require(GEOquery)
require(org.Hs.eg.db)

rootDir <- "/home/shraddhapai/BaderLab/2018_AsthmaPBMC"
inDir <- sprintf("%s/input",rootDir)
datDir <- sprintf("%s/output/basic_180220/pred",rootDir)
outDir <- sprintf("%s/output/basic_180220/plot",rootDir)
pathFile <-sprintf("%s/anno/Human_AllPathways_November_01_2017_symbol.gmt",
	rootDir)

dat <- getGEO(filename=sprintf("%s/GSE40732_series_matrix.txt.gz",inDir),
	GSEMatrix=TRUE)
xpr <- exprs(dat)
pheno <- pData(dat)
# map GB ID to symbol
x <- mapIds(org.Hs.eg.db, keys=rownames(xpr), 
	column="SYMBOL",keytype="ACCNUM",
	multiVals="first")
common <- intersect(names(x),rownames(xpr))
xpr <- xpr[which(rownames(xpr) %in% common),]
x <- x[which(names(x) %in% common)]

midx <- match(rownames(xpr),names(x))
gnames <- x[midx]
agg <- aggregate(xpr, by=list(gene_name=gnames),FUN=mean)
xpr <- agg[,-1]
rownames(xpr) <- agg[,1]

pheno <- pheno[,c("geo_accession","characteristics_ch1")]
st <- rep(NA, nrow(pheno))
st[which(pheno[,2] %in% "asthma: FALSE")] <- "control"
st[which(pheno[,2] %in% "asthma: TRUE")] <- "asthma"
pheno[,2] <- st
colnames(pheno) <- c("ID","STATUS")
pheno[,1] <- as.character(pheno[,1])

pathwayList <- readPathways(pathFile)
head(pathwayList)

if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
#predDir <- sprintf("%s/rng%i/predictionResults.txt", 
#	datDir,1:11)
predPerf <- plotPerf(datDir, predClasses=predClasses)
dev.off()
auroc <- unlist(lapply(predPerf, function(x) x$auroc))
aupr <- unlist(lapply(predPerf, function(x) x$aupr))
acc <- unlist(lapply(predPerf, function(x) x$accuracy))

###featList <- list(
###	control=sprintf("%s/rng%i/control/GM_results/control_pathway_CV_score.txt",
###		datDir,1:11),
###	asthma=sprintf("%s/rng%i/asthma/GM_results/asthma_pathway_CV_score.txt",
###		datDir,1:11))
featScores <- getFeatureScores(datDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.9)
})

netInfoFile <- sprintf("%s/inputNets.txt",datDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir)
