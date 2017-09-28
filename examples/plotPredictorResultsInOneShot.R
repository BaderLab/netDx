## ------------------------------------------------------------------------
#httr
tryCatch(expr = { library(httr)}, 
          error = function(e) { install.packages("httr")}, finally = library(httr))

#RJSONIO
tryCatch(expr = { library(RJSONIO)}, 
          error = function(e) { install.packages("RJSONIO")}, finally = library(RJSONIO))

#r2cytoscape
tryCatch(expr = { library(r2cytoscape)}, 
          error = function(e) { devtools::install_github('cytoscape/cytoscape-automation/for-scripters/R/r2cytoscape')}, finally = library(r2cytoscape))

# EasycyRest
tryCatch(expr = { library(EasycyRest); detach(package:EascyRest,unload=TRUE)}, 
          error = function(e) { devtools::install_github('BaderLab/Easycyrest/EasycyRest@0.1')}, finally = {})

## ----eval=TRUE-----------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))
suppressWarnings(suppressMessages(require(netDx.examples)))

## ----eval=TRUE-----------------------------------------------------------
phenoFile <- sprintf("%s/extdata/KIRC_pheno.rda",path.package("netDx.examples"))
lnames <- load(phenoFile)
head(pheno)

## ------------------------------------------------------------------------
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

## ------------------------------------------------------------------------
xpr_genes <- sprintf("%s/extdata/EMap_input/genenames.txt",
      path.package("netDx.examples"))
xpr_genes <- read.delim(xpr_genes,h=FALSE,as.is=TRUE)[,1]
head(xpr_genes)

## ------------------------------------------------------------------------
pathwayList <- lapply(pathwayList, function(x) x[which(x %in% xpr_genes)])

## ------------------------------------------------------------------------
inDir <- sprintf("%s/extdata/KIRC_output",
	path.package("netDx.examples"))
out <- plotAllResults(pheno, inDir,outDir=sprintf("%s/plots",getwd()),
               fsCutoff=10,fsPctPass=0.7,pathwaySet=pathwayList)

## ------------------------------------------------------------------------
dir(sprintf("%s/plots",getwd()),recursive=TRUE)

## ------------------------------------------------------------------------
sessionInfo()

