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

## ------------------------------------------------------------------------
# Basic settings
version.url <-"http://localhost:1234/v1/version"
cytoscape.open <- TRUE

tryCatch(expr = { GET(version.url)}, 
         error = function(e) { return (cytoscape.open = FALSE)}, finally =function(r){ return(cytoscape.open = TRUE)})
         
if(!cytoscape.open){
  #try and launch cytoscape
 stop("Cytoscape is not open.  Please launch cytoscape.")
} else{
  cytoscape.version <-  GET(version.url)
  cy.version <- fromJSON(rawToChar(cytoscape.version$content))
  print(cy.version)
}

## ----eval=FALSE----------------------------------------------------------
## outerLoop <- 100     # num times to split data into train/blind test samples
## innerLoop <- 10      # num folds for cross-validation, also max score for a network
## netScores <- list()  # collect <outerLoop> set of netScores
## perf <- list()       # collect <outerLoop> set of test evaluations
## 
## for k in 1:outerLoop
##  [train, test] <- splitData(80:20) # split data using RNG seed
##  netScores[[k]] <- runCV(train)
##  perf[[k]] <- collectPerformance(netScores[[k]], test)
## end

## ----eval=TRUE-----------------------------------------------------------
suppressMessages(require(netDx))
suppressMessages(require(netDx.examples))

## ----eval=TRUE-----------------------------------------------------------
phenoFile <- sprintf("%s/extdata/KIRC_pheno.rda",path.package("netDx.examples"))
lnames <- load(phenoFile)
head(pheno)

## ------------------------------------------------------------------------
outDir <- paste(getwd(),"plots",sep="/")
if (!file.exists(outDir)) dir.create(outDir)
setwd(outDir)

## ------------------------------------------------------------------------
inDir <- sprintf("%s/extdata/KIRC_output",
	path.package("netDx.examples"))
all_rngs <- list.dirs(inDir, recursive = FALSE)
print(head(basename(all_rngs)))

## ----eval=TRUE-----------------------------------------------------------
dir(all_rngs[1])

## ----eval=TRUE,fig.width=6,fig.height=6----------------------------------
predClasses <- c("SURVIVEYES","SURVIVENO")
predFiles <- unlist(lapply(all_rngs, function(x) 
		paste(x, "predictionResults.txt", sep = "/")))
pdf(sprintf("%s/perf.pdf",outDir),height=8,width=6)
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()

## ----eval=TRUE-----------------------------------------------------------
featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))

## ----eval=TRUE-----------------------------------------------------------
dim(featScores[[1]])

## ----eval=TRUE-----------------------------------------------------------
head(featScores[[1]][,1:10])

## ----eval=TRUE-----------------------------------------------------------
featSelNet <- lapply(featScores, function(x) {
	callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
})

## ----eval=TRUE-----------------------------------------------------------
tmp <- lapply(featSelNet,print)

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
netInfoFile <- sprintf("%s/extdata/KIRC_output/inputNets.txt",
      path.package("netDx.examples"))
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
head(netInfo)

## ------------------------------------------------------------------------
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir)

## ------------------------------------------------------------------------
pngFiles <- list()
for (curGroup in names(EMap_input)[1:2]) {
	pngFiles[[curGroup]] <- plotEmap(gmtFile=EMap_input[[curGroup]][1], 
		                        nodeAttrFile=EMap_input[[curGroup]][2],
		                        netName=curGroup,outDir=outDir)
}

## ----eval=TRUE-----------------------------------------------------------
netInfo <- plotIntegratedPSN(pheno=pheno,baseDir=sprintf("%s/rng1",inDir),
	netNames=featSelNet,outDir=outDir)

## ----eval=TRUE-----------------------------------------------------------
summary(netInfo)

## ------------------------------------------------------------------------
sessionInfo()

