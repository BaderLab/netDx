rm(list=ls())
suppressWarnings(suppressMessages(require(netDx)))

library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
curatedTCGAData(diseaseCode="BRCA", assays="*",dru.run=TRUE)

# fetch mrna, clinical data
brca <- curatedTCGAData("BRCA",c("mRNAArray"),FALSE)

# prepare data 
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

idx <- union(which(pam50 == "Normal-like"), which(is.na(st2)))
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

# remove duplicate arrays
smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- brca[[1]][,notdup]

# create "ID" and "STATUS" column netDx looks for
# These should be columns of the sample metadata in colData()
colData(brca)$STATUS <- pam50
pID <- colData(brca)$patientID
colData(brca)$ID <- pID


groupList <- list()
pathList <- readPathways(getExamplePathways())
groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
groupList[["clinical"]] <- list(age="patient.age_at_initial_pathologic_diagnosis",
	stage="STAGE")

## ------------------------------------------------------------------------
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) { 
	netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
				rownames(dataList[["BRCA_mRNAArray-20160128"]]),
			   	groupList[["BRCA_mRNAArray-20160128"]],
				netDir,verbose=FALSE, 
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
		simMetric="custom",customFunc=normDiff, # custom function
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
out <- buildPredictor(dataList=brca,groupList=groupList,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred_output",tempdir()), ## absolute path
   numCores=1L,featScoreMax=2L, featSelCutoff=1L,numSplits=2L)

print(summary(out))



## ------------------------------------------------------------------------
sessionInfo()

