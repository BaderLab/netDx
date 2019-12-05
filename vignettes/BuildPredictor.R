## ----eval=FALSE----------------------------------------------------------
## 
## # load libraries
## suppressWarnings(suppressMessages(require(netDx)))
## options(stringsAsFactors = FALSE)
## 
## # prepare data
## library(curatedTCGAData)
## library(MultiAssayExperiment)
## curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)
## brca <- curatedTCGAData("BRCA",c("mRNAArray"),FALSE)
## 
## staget <- colData(brca)$pathology_T_stage
## st2 <- rep(NA,length(staget))
## st2[which(staget %in% c("t1","t1a","t1b","t1c"))] <- 1
## st2[which(staget %in% c("t2","t2a","t2b"))] <- 2
## st2[which(staget %in% c("t3","t3a"))] <- 3
## st2[which(staget %in% c("t4","t4b","t4d"))] <- 4
## colData(brca)$STAGE <- st2
## 
## pam50 <- colData(brca)$PAM50.mRNA
## pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
## pam50[which(pam50 %in% "Luminal A")] <- "LumA"
## colData(brca)$pam_mod <- pam50
## 
## idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
## pID <- colData(brca)$patientID
## tokeep <- setdiff(pID, pID[idx])
## brca <- brca[,tokeep,]
## 
## smp <- sampleMap(brca)
## samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
## notdup <- samps[which(!duplicated(samps$primary)),"colname"]
## brca[[1]] <- brca[[1]][,notdup]
## 
## # set "ID" and "STATUS" columns (netDx looks for these).
## pID <- colData(brca)$patientID
## colData(brca)$ID <- pID
## colData(brca)$STATUS <- colData(brca)$pam_mod
## 
## # define features
## groupList <- list()
## 
## # genes in mRNA data are grouped by pathways
## pathList <- readPathways(getExamplePathways())
## groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
## # clinical data is not grouped; each variable is its own feature
## groupList[["clinical"]] <- list(
##       age="patient.age_at_initial_pathologic_diagnosis",
## 	   stage="STAGE"
## )
## 
## # define simliarity function used to create features
## # in this example, pairwise Pearson correlation is used for gene expression
## # and normalized difference is used for clinical data
## makeNets <- function(dataList, groupList, netDir,...) {
## 	netList <- c() # initialize before is.null() check
## 	# make RNA nets (NOTE: the check for is.null() is important!)
## 	# (Pearson correlation)
## 	if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) {
## 	netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
## 				rownames(dataList[["BRCA_mRNAArray-20160128"]]),
## 			   	groupList[["BRCA_mRNAArray-20160128"]],
## 				netDir,verbose=FALSE,
## 			  	writeProfiles=TRUE,...)
## 	}
## 	
## 	# make clinical nets (normalized difference)
## 	netList2 <- c()
## 	if (!is.null(groupList[["clinical"]])) {
## 	netList2 <- makePSN_NamedMatrix(dataList$clinical,
## 		rownames(dataList$clinical),
## 		groupList[["clinical"]],netDir,
## 		simMetric="custom",customFunc=normDiff, # custom function
## 		writeProfiles=FALSE,
## 		sparsify=TRUE,verbose=TRUE,append=TRUE,...)
## 	}
## 	netList <- c(unlist(netList),unlist(netList2))
## 	return(netList)
## }
## 
## # train the model.
## # Here we run two train/test splits (numSplits). In each split,
## # feature selection scores features out of 2, and features that
## # score >=1 are used to classify test samples
## 
## set.seed(42) # make results reproducible
## out <- buildPredictor(dataList=brca,groupList=groupList,
##    makeNetFunc=makeNets, ### custom network creation function
##    outDir=sprintf("%s/pred_output",tempdir()), ## absolute path
##    numCores=1L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)
## 
## # look at results
## print(summary(out))
## 


## ----eval=FALSE----------------------------------------------------------
## numSplits <- 100     # num times to split data into train/blind test samples
## featScoreMax <- 10      # num folds for cross-validation, also max score for a network
## featSelCutoff <- 9
## netScores <- list()  # collect <numSplits> set of netScores
## perf <- list()       # collect <numSplits> set of test evaluations
## 
## for k in 1:numSplits
##  [train, test] <- splitData(80:20) # split data using RNG seed
##   featScores[[k]] <- scoreFeatures(train, featScoreMax)
##  topFeat[[k]] <- applyFeatCutoff(featScores[[k]])
##  perf[[k]] <- collectPerformance(topFeat[[k]], test)
## end


## ------------------------------------------------------------------------
suppressWarnings(suppressMessages(require(netDx)))


## ----eval=TRUE-----------------------------------------------------------
suppressMessages(library(curatedTCGAData))
suppressMessages(library(MultiAssayExperiment))


## ------------------------------------------------------------------------
curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)


## ------------------------------------------------------------------------
brca <- curatedTCGAData("BRCA",c("mRNAArray"),FALSE)


## ----eval=TRUE-----------------------------------------------------------
staget <- colData(brca)$pathology_T_stage
st2 <- rep(NA,length(staget))
st2[which(staget %in% c("t1","t1a","t1b","t1c"))] <- 1
st2[which(staget %in% c("t2","t2a","t2b"))] <- 2
st2[which(staget %in% c("t3","t3a"))] <- 3
st2[which(staget %in% c("t4","t4b","t4d"))] <- 4
colData(brca)$STAGE <- st2

pam50 <- colData(brca)$PAM50.mRNA
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
colData(brca)$pam_mod <- pam50

idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
pID <- colData(brca)$patientID
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- brca[[1]][,notdup]


## ----eval=TRUE-----------------------------------------------------------
pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- colData(brca)$pam_mod


## ----eval=TRUE-----------------------------------------------------------
ids <- sprintf("patient%i",1:20)
mrna <- matrix(rnorm(2000),nrow=100,ncol=20) # 100 genes x 20 patients
rownames(mrna) <- sprintf("gene%i",1:100)
colnames(mrna) <- ids

age <- round(runif(20,min=20,max=35))
important_variable <- c(rep("LOW",10),rep("HIGH",10))
clin <- t(data.frame(age=age,imp_var=important_variable))
colnames(clin) <- ids

dataList <- list(clinical=clin,transcription=mrna)

summary(dataList)


## ----eval=TRUE-----------------------------------------------------------
groupList <- list()

# genes in mRNA data are grouped by pathways
pathList <- readPathways(getExamplePathways())
groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
# clinical data is not grouped; each variable is its own feature
groupList[["clinical"]] <- list(
      age="patient.age_at_initial_pathologic_diagnosis",
	   stage="STAGE"
)


## ----eval=TRUE-----------------------------------------------------------
summary(groupList)


## ----eval=TRUE-----------------------------------------------------------
groupList[["BRCA_mRNAArray-20160128"]][1:3]


## ----eval=TRUE-----------------------------------------------------------
head(groupList[["clinical"]])


## ------------------------------------------------------------------------
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
		sparsify=TRUE,verbose=TRUE,...)
	}
	netList <- c(unlist(netList),unlist(netList2))
	return(netList)
}



## ----eval=TRUE-----------------------------------------------------------
set.seed(42) # make results reproducible
outDir <- sprintf("%s/pred_output",tempdir()) # location for intermediate work
# set keepAllData to TRUE to not delete at the end of the predictor run.
# This can be useful for debugging.

out <- buildPredictor(dataList=brca,groupList=groupList,
  makeNetFunc=makeNets,outDir=outDir,
  numSplits=2L,featScoreMax=2L, featSelCutoff=1L,
  keepAllData=2,
	numCores=1L)


## ----eval=TRUE-----------------------------------------------------------
summary(out)
summary(out$Split1)


## ----eval=TRUE-----------------------------------------------------------
save(out,file=sprintf("%s/results.rda",outDir))


## ----eval=TRUE-----------------------------------------------------------
numSplits <- 2
st <- unique(colData(brca)$STATUS) # to get similarity scores for each class
for (k in 1:numSplits) { 
	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
	oF <- sprintf("%s/Split%i_predictionResults.txt",outDir,k)
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",sprintf("%s_SCORE",st))]
	write.table(tmp,file=oF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)
}


## ------------------------------------------------------------------------
sessionInfo()

