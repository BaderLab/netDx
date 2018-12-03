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

## ------------------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))
suppressWarnings(suppressMessages(require(netDx.examples)))

## ----eval=T--------------------------------------------------------------
load(sprintf("%s/extdata/nestedCV_input.rda",
             path.package("netDx.examples")))
head(pheno)

## ----eval=TRUE-----------------------------------------------------------
names(dats)

## ----eval=TRUE-----------------------------------------------------------
head(dats[["rna"]][,1:6])

## ----eval=TRUE-----------------------------------------------------------
head(dats[["clinical"]][,1:6])

## ----eval=TRUE-----------------------------------------------------------
lapply(dats, nrow)

## ----eval=TRUE-----------------------------------------------------------
names(groupList)

## ----eval=TRUE-----------------------------------------------------------
groupList[["rna"]][1:3]

## ----eval=TRUE-----------------------------------------------------------
head(groupList[["clinical"]])

## ----eval=TRUE-----------------------------------------------------------
groupList[["rna"]] <- groupList[["rna"]][1:3]

## ------------------------------------------------------------------------
KIRC_makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["rna"]])) { 
	netList <- makePSN_NamedMatrix(dataList$rna, 
					rownames(dataList$rna),
			   	groupList[["rna"]],netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	netList <- unlist(netList)
	cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
	}
	
	# make clinical nets,one net for each variable
	netList2 <- c()
	if (!is.null(groupList[["clinical"]])) {
	netList2 <- makePSN_NamedMatrix(dataList$clinical, 
		rownames(dataList$clinical),
		groupList[["clinical"]],netDir,
		simMetric="custom",customFunc=netDx::normDiff, # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,verbose=TRUE,append=TRUE,...)
	}
	netList2 <- unlist(netList2)
	cat(sprintf("Made %i clinical nets\n", length(netList2)))
	netList <- c(netList,netList2) 
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}

## ----eval=TRUE-----------------------------------------------------------
runPredictor_nestedCV(pheno,
   dataList=dats,groupList=groupList,
   makeNetFunc=KIRC_makeNets, ### custom network creation function
   outDir=sprintf("%s/nestedCV_output",getwd()), ## absolute path
   numCores=1L,nFoldCV=2L, CVcutoff=1L,numSplits=2L)

## ----eval=TRUE-----------------------------------------------------------
outDir <- sprintf("%s/nestedCV_output",getwd())
dir(outDir)

## ------------------------------------------------------------------------
dir(sprintf("%s/rng1",outDir))
pred <- read.delim(sprintf("%s/rng1/predictionResults.txt",outDir),h=T,as.is=T)
head(pred)

## ------------------------------------------------------------------------
dir(sprintf("%s/rng1/SURVIVEYES",outDir))
sc <- read.delim(sprintf("%s/rng1/SURVIVEYES/GM_results/SURVIVEYES_pathway_CV_score.txt",outDir),
   sep="\t",h=T,as.is=T)
head(sc)

## ------------------------------------------------------------------------
sessionInfo()

