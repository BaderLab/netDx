## ----eval=TRUE-----------------------------------------------------------

# load libraries
suppressWarnings(suppressMessages(require(netDx)))

# prepare data
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)
brca <- curatedTCGAData("BRCA",c("mRNAArray"),FALSE)

pam50 <- colData(brca)$PAM50.mRNA
staget <- colData(brca)$pathology_T_stage
st2 <- rep(NA,length(staget))
st2[which(staget %in% c("t1","t1a","t1b","t1c"))] <- 1
st2[which(staget %in% c("t2","t2a","t2b"))] <- 2
st2[which(staget %in% c("t3","t3a"))] <- 3
st2[which(staget %in% c("t4","t4b","t4d"))] <- 4
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
colData(brca)$STAGE <- st2
colData(brca)$pam_mod <- pam50

## ----eval=TRUE-----------------------------------------------------------
idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- brca[[1]][,notdup]

# set "ID" and "STATUS" columns (netDx looks for these). 
pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- colData(brca)$pam_mod


## ----eval=TRUE-----------------------------------------------------------
# define features 
groupList <- list()

# genes in mRNA data are grouped by pathways
pathList <- readPathways(getExamplePathways())
groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
# clinical data is not grouped; each variable is its own feature
groupList[["clinical"]] <- list(
      age="patient.age_at_initial_pathologic_diagnosis",
	   stage="STAGE"
)

# define simliarity function used to create features
# in this example, pairwise Pearson correlation is used for gene expression
# and normalized difference is used for clinical data
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() # initialize before is.null() check
	# make RNA nets (NOTE: the check for is.null() is important!)
	# (Pearson correlation)
	if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) { 
	netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
				rownames(dataList[["BRCA_mRNAArray-20160128"]]),
			   	groupList[["BRCA_mRNAArray-20160128"]],
				netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	}
	
	# make clinical nets (normalized difference)
	netList2 <- c()
	if (!is.null(groupList[["clinical"]])) {
	netList2 <- makePSN_NamedMatrix(dataList$clinical, 
		rownames(dataList$clinical),
		groupList[["clinical"]],netDir,
		simMetric="custom",customFunc=normDiff, # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,verbose=TRUE,append=TRUE,...)
	}
	netList <- c(unlist(netList),unlist(netList2))
	return(netList)
}

# train the model. 
# Here we run two train/test splits (numSplits). In each split, 
# feature selection scores features out of 2, and features that
# score >=1 are used to classify test samples

set.seed(42) # make results reproducible
out <- buildPredictor(dataList=brca,groupList=groupList,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred_output",tempdir()), ## absolute path
   numCores=1L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)

# look at results
print(summary(out))



## ------------------------------------------------------------------------
sessionInfo()



## ----eval=FALSE----------------------------------------------------------
## numSplits <- 100     # num times to split data into train/blind test samples
## featScoreMax <- 10      # num folds for cross-validation, also max score for a network
## netScores <- list()  # collect <numSplits> set of netScores
## perf <- list()       # collect <numSplits> set of test evaluations
## 
## for k in 1:numSplits
##  [train, test] <- splitData(80:20) # split data using RNG seed
##  netScores[[k]] <- scoreFeatures(train, featScoreMax)
##  perf[[k]] <- collectPerformance(netScores[[k]], test)
## end


## ------------------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))
suppressWarnings(suppressMessages(require(netDx.examples)))


## ----eval=FALSE----------------------------------------------------------
## load(sprintf("%s/extdata/buildPred_input.rda",
##              path.package("netDx.examples")))
## head(pheno)


## ----eval=FALSE----------------------------------------------------------
## names(dats)


## ----eval=FALSE----------------------------------------------------------
## head(dats[["rna"]][,1:6])


## ----eval=FALSE----------------------------------------------------------
## head(dats[["clinical"]][,1:6])


## ----eval=FALSE----------------------------------------------------------
## lapply(dats, nrow)


## ----eval=FALSE----------------------------------------------------------
## names(groupList)


## ----eval=FALSE----------------------------------------------------------
## groupList[["rna"]][1:3]


## ----eval=FALSE----------------------------------------------------------
## head(groupList[["clinical"]])


## ----eval=FALSE----------------------------------------------------------
## groupList[["rna"]] <- groupList[["rna"]][1:3]


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


## ----eval=FALSE----------------------------------------------------------
## buildPredictor(pheno,
##    dataList=dats,groupList=groupList,
##    makeNetFunc=KIRC_makeNets, ### custom network creation function
##    outDir=sprintf("%s/pred_output",tempdir()), ## absolute path
##    numCores=1L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)


## ----eval=FALSE----------------------------------------------------------
## outDir <- sprintf("%s/pred_output",tempdir())
## dir(outDir)


## ----eval=FALSE----------------------------------------------------------
## dir(sprintf("%s/rng1",outDir))
## pred <- read.delim(sprintf("%s/rng1/predictionResults.txt",outDir),h=TRUE,as.is=TRUE)
## head(pred)


## ------------------------------------------------------------------------
sessionInfo()

