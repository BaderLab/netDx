#' feature selection for GBM from PanCancer survival dataset
#' 10-fold CV predictor design with clinical and mRNA data
rm(list=ls())
require(netDx)
require(netDx.examples)

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.8
cutoff <- 9

phenoFile <- sprintf("%s/extdata/KIRC_pheno.rda",
	path.package("netDx.examples"))
lnames <- load(phenoFile)
head(pheno)

# filter pathwayList for EMap
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
xpr_genes <- sprintf("%s/extdata/EMap_input/genenames.txt",
      path.package("netDx.examples"))
xpr_genes <- read.delim(xpr_genes,h=FALSE,as.is=TRUE)[,1]

inDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/pathway_170502"
outDir <- sprintf("%s/results",inDir)

res <- plotAllResults(pheno=pheno,inDir=inDir,outDir,pathwaySet=pathwayList)
