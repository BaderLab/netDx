# testing BRCA pathway-level correlation view
rm(list=ls())

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

subtypes<- c("LumA")
pheno$STATUS[which(!pheno$STATUS %in% subtypes)] <- "other"
pheno$STATUS <- factor(pheno$STATUS, levels=c("other","LumA"))
subtypes <- c(subtypes,"other") # add residual

rootDir <- "/home/shraddhapai/BaderLab/2017_BRCA"
pathFile <- sprintf("%s/anno/Human_AllPathways_February_01_2018_symbol.gmt",
	rootDir)
pathwayList <- readPathways(pathFile)
head(pathwayList)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

inputNets <- data.frame(datatype="rna",netName=names(pathwayList))

inDir <- sprintf("%s/output/BRCA_part2_180223",rootDir)
predClasses <- c("LumA","other")
featScores <- getFeatureScores(inDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
})

###selFeat <- c(
###	"ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS.profile",
###	"ADENOSINE_DEOXYRIBONUCLEOTIDES__I_DE_NOVO__I__BIOSYNTHESIS.profile",
###	"MISMATCH_REPAIR__MMR__DIRECTED_BY_MSH2:MSH3__MUTSBETA_.profile",
###	"MITOTIC_SPINDLE_CHECKPOINT.profile",
###	"RESOLUTION_OF_D-LOOP_STRUCTURES.profile"
###)
set.seed(123)
selFeat <- featSelNet[["LumA"]]
selFeat <- selFeat[sample(1:length(selFeat),10,replace=F)]
selFeat  <- c(selFeat, "MISMATCH_REPAIR__MMR__DIRECTED_BY_MSH2:MSH3__MUTSBETA_.profile")
selFeat <- sub(".profile","",selFeat)

source("multiplot.R")
require(plotrix)
require(ggplot2)
source("corrFeatWithOutcome.R");
resMat <- corrFeatWithOutcome(pheno,dats,gps,inputNets,selFeat,numPCs=3,filePfx="LumA_top10",plotLevels=c("LumA","other"))

