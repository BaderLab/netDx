#' plot BRCA results

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

rootDir <- "/home/shraddhapai/BaderLab/2017_BRCA"

pathFile <- sprintf("%s/anno/Human_AllPathways_February_01_2018_symbol.gmt",
	rootDir)
pathwayList <- readPathways(pathFile)
head(pathwayList)

xpr_genes <- rownames(xpr)
pathwayList <- lapply(pathwayList,function(x) x[which(x %in% xpr_genes)])

pheno$STATUS[which(pheno$STATUS!="LumA")] <- "other"

#out <- plotAllResults(pheno, sprintf("%s/pred",rootDir),
#							 outDir=sprintf("%s/plot",rootDir),
#               fsCutoff=10,fsPctPass=0.7,pathwaySet=pathwayList)

inDir <- sprintf("%s/output/BRCA_part2_180223",rootDir)
outDir <- sprintf("%s/output/BRCA_part2_180223/plot",rootDir)
if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))

#predFiles <- sprintf("%s/rng%i/predictionResults.txt", inDir,1:79)
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()

featFiles <- list(
	LumA=sprintf("%s/rng%i/LumA/GM_results/LumA_pathway_CV_score.txt", inDir,1:79),
	other=sprintf("%s/rng%i/other/GM_results/other_pathway_CV_score.txt", inDir,1:79)
)
#featScores <- getFeatureScores(featFiles,predClasses=predClasses)
featScores <- getFeatureScores(inDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
})

netInfoFile <- sprintf("%s/inputNets.txt",inDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir,pctPass=0.7)

auroc <- unlist(lapply(predPerf,function(x) x$auroc))
aupr <- unlist(lapply(predPerf,function(x) x$aupr))
acc <- unlist(lapply(predPerf,function(x) x$accuracy))
cat(sprintf("mean auroc = %1.2f ; aupr = %1.2f ; acc = %1.2f%%",
	mean(auroc), mean(aupr), mean(acc)))

cat("Performance\n")
#cat(sprintf("AUROC = %1.2f +/- %1.2f\n", mean(auroc),sd(auroc)/sqrt(length(auroc))))
#cat(sprintf("AUPR = %1.2f +/- %1.2f\n", mean(aupr),sd(aupr)/sqrt(length(aupr))))
#cat(sprintf("ACC = %1.2f +/- %1.2f\n", mean(acc),sd(acc)/sqrt(length(acc))))
cat(sprintf("AUROC = %1.2f +/- %1.2f\n", mean(auroc),sd(auroc)))
cat(sprintf("AUPR = %1.2f +/- %1.2f\n", mean(aupr),sd(aupr)))
cat(sprintf("ACC = %1.2f +/- %1.2f\n", mean(acc),sd(acc)))
