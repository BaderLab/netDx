## ----eval=TRUE----------------------------------------------------------------


## ---- fig.cap="Lab 1 design: We will integrate four layers of genomic data. Each layer will be converted into a single patient similarity network using Pearson correlation for pairwise similarity. ",echo=FALSE----
knitr::include_graphics("images/vignette1_design.jpg")


## ---- class.source="codeblock",eval=TRUE--------------------------------------
suppressMessages(library(curatedTCGAData))


## ----eval=TRUE----------------------------------------------------------------
suppressMessages(library(curatedTCGAData))


## ---- class.source="codeblock",eval=TRUE--------------------------------------
curatedTCGAData(diseaseCode="BRCA", assays="*",dry.run=TRUE, version="1.1.38")


## ---- class.source="codeblock",eval=TRUE--------------------------------------
brca <- suppressMessages(curatedTCGAData("BRCA",
                                         c("mRNAArray","Methylation_methyl27", 
										 "RPPAArray","miRNASeqGene"),
                                         dry.run=FALSE, version="1.1.38"))


## ---- class.source="codeblock",eval=TRUE--------------------------------------
brca


## ---- class.source="codeblock",eval=TRUE--------------------------------------
summary(assays(brca))


## ---- class.source="codeblock",eval=TRUE--------------------------------------
names(assays(brca))


## ---- class.source="codeblock",eval=TRUE--------------------------------------
mir <- assays(brca)[["BRCA_miRNASeqGene-20160128"]]
head(mir[,1:5])


## ---- class.source="codeblock",eval=TRUE--------------------------------------
pheno <- colData(brca)
colnames(pheno)[1:20]
head(pheno[,1:5])


## ---- class.source="codeblock",eval=TRUE--------------------------------------
source("prepare_data.R")
brca <- prepareData(brca)


## ---- class.source="codeblock",eval=TRUE--------------------------------------
pheno <- colData(brca)
head(pheno[,c("ID","STATUS")])
table(pheno$STATUS,useNA="always")  # good practice: useNA="always" shows missing values


## ---- class.source="codeblock",eval=TRUE--------------------------------------
suppressWarnings(suppressMessages(require(netDx)))


## ---- class.source="codeblock",eval=TRUE--------------------------------------
expr <- assays(brca)
groupList <- list()
for (k in 1:length(expr)) {	# loop over all layers
	cur <- expr[[k]]; nm <- names(expr)[k]

	# all measures from this layer go into our single PSN
	groupList[[nm]] <- list(nm=rownames(cur)) 

	# assign same layer name as in input data
	names(groupList[[nm]])[1] <- nm;
}


## ---- class.source="codeblock",eval=TRUE--------------------------------------
summary(groupList)


## ---- class.source="codeblock",eval=TRUE--------------------------------------
names(groupList[["BRCA_mRNAArray-20160128"]])
length(groupList[["BRCA_mRNAArray-20160128"]][[1]])
head(groupList[["BRCA_mRNAArray-20160128"]][[1]])


## ---- class.source="codeblock", eval=TRUE-------------------------------------
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() # initialize before is.null() check
	
	layerNames <- c("BRCA_miRNASeqGene-20160128",
		"BRCA_mRNAArray-20160128",
		"BRCA_RPPAArray-20160128",
		"BRCA_Methylation_methyl27-20160128")
	
	for (nm in layerNames){  			## for each layer
		if (!is.null(groupList[[nm]])){ ## must check for null for each layer
			netList_cur <- makePSN_NamedMatrix(
				dataList[[nm]],
				rownames(dataList[[nm]]),	## names of measures (e.g. genes, CpGs)
				groupList[[nm]],			## how to group measures in that layer
				netDir,						## leave this as-is, netDx will figure out where this is.
				verbose=FALSE, 			
				writeProfiles=TRUE,   		## use Pearson correlation-based similarity
				...
				)

			netList <- c(netList,netList_cur)	## just leave this in
		}
	}
	return(unlist(netList))	## just leave this in 
}


## In practice a good starting point is `featScoreMax=10`, `featSelCutoff=9` and `numSplits=10L`, but these parameters depend on the sample sizes in the dataset and heterogeneity of the samples.


## ----lab1-buildpredictor ,class.source="codeblock",eval=TRUE------------------
t0 <- Sys.time()
set.seed(42) # make results reproducible
outDir <- paste(tempdir(),randAlphanumString(),
	"pred_output",sep=getFileSep())
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
model <- suppressMessages(buildPredictor(
	dataList=brca,			## your data
	groupList=groupList,	## grouping strategy
	makeNetFunc=makeNets,	## function to build PSNs
	outDir=outDir, 			## output directory
	trainProp=0.8,			## pct of samples to use to train model in
							## each split
	numSplits=2L,			## number of train/test splits
	featSelCutoff=1L,		## threshold for calling something
							## feature-selected
	featScoreMax=2L,		## max score for feature selection
 numCores=2L,			## set higher for parallelizing
 debugMode=FALSE,
 keepAllData=FALSE,	## set to TRUE for debugging or low-level files used by the dictor
 logging="none"
  ))
t1 <- Sys.time()
print(t1-t0)

browser()

## ---- class.source="codeblock",eval=TRUE--------------------------------------
outFile <- sprintf("%s/CBW_Lab1_full.rda",tempdir())
###download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/Lab1_files/Lab1_10splits.rda",
###	destfile=outFile)
lnames <- load(outFile)


## ----lab1-getresults,class.source="codeblock",eval=TRUE-----------------------
results <- getResults(brca,model_full,featureSelCutoff=9L,
	featureSelPct=0.9)


## ---- class.source="codeblock",eval=TRUE--------------------------------------
summary(results)


## ---- class.source="codeblock",eval=TRUE--------------------------------------
results$performance


## ---- class.source="codeblock", eval=TRUE-------------------------------------
results$featureScores


## ---- class.source="codeblock",eval=TRUE--------------------------------------
confMat <- confusionMatrix(model_full)


## ---- class.source="codeblock",eval=TRUE--------------------------------------
results$selectedFeatures


## ---- class.source="codeblock",fig.width=8,fig.height=8, eval=TRUE------------
## this call doesn't work in Rstudio; for now we've commented this out and saved the PSN file. 
## psn <- getPSN(brca,groupList_full,makeNets_full,results$selectedFeatures)

psnFile <- sprintf("%s/psn.rda",tempdir())
download.file("https://github.com/RealPaiLab/CBW_CAN_DataIntegration_2021/raw/master/supporting_files/Lab1_files/Lab1_PSN.rda",
	destfile=outFile)
load(outFile)

require(Rtsne)
tsne <- tSNEPlotter(
	psn$patientSimNetwork_unpruned, 
	colData(brca)
	)


## -----------------------------------------------------------------------------
sessionInfo()

