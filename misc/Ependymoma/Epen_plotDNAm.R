#' plot BRCA results
rm(list=ls())
require(netDx)
require(netDx.examples)
#rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_Ependymoma"
rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"

phenoFile <- "/home/shraddhapai/BaderLab/2018_Epen_DNAm/input/GSE90496_pData.txt"
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
ttype <- pheno$characteristics_ch1

inFile <- sprintf("%s/input/netDx_prepared/Ependymoma_cohortMerged_180125.Rdata",rootDir)
load(inFile)
# ----------------------
# input processing
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
ttype <- pheno$characteristics_ch1
idx <- which(ttype %in% c("methylation class: EPN, PF A","methylation class: EPN, PF B"))
cat(sprintf("Got %i samples\n",length(idx)))
pheno <- pheno[idx,] # limit to EPN samples
cpos <- regexpr("sample", pheno$title)
bpos <- regexpr("\\[reference", pheno$title)
str <- as.integer(substr(pheno$title, cpos+7, bpos-2)) # get sample number
pheno$ID <- paste("SAMPLE", str,sep=".")
pheno <- pheno[,c("ID","characteristics_ch1")]
st <- rep("",nrow(pheno))
st[grep("PF A", pheno[,2])] <- "PFA"
st[grep("PF B", pheno[,2])] <- "PFB"
pheno$STATUS <- st

#out <- plotAllResults(pheno, sprintf("%s/pred",rootDir),
#							 outDir=sprintf("%s/plot",rootDir),
#               fsCutoff=10,fsPctPass=0.7,pathwaySet=pathwayList)

#setName <- "Epen_prunedOneNet_0.001_180409"
#setName <- "Epen_prunedPathway_0.01_180409"
setName <- "Epen_prunedOneNet_0.001_180410"
inDir <- sprintf("%s/output/%s/pred",rootDir,setName)
outDir <- sprintf("%s/output/%s/plot",rootDir,setName)

if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()
auroc <- unlist(lapply(predPerf, function(x) x$auroc*100))
aupr <- unlist(lapply(predPerf, function(x) x$aupr*100))
acc <- unlist(lapply(predPerf, function(x) x$accuracy))

cat("--------------\n")
cat(sprintf("Performance: %s\n",setName))
cat(sprintf("AUROC = %1.2f +/- %1.2f\n",mean(auroc),sd(auroc)))
cat(sprintf("AUPR = %1.2f +/- %1.2f\n",mean(aupr),sd(aupr)))
cat(sprintf("Accuracy = %1.2f +/- %1.2f\n",mean(acc),sd(acc)))
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
