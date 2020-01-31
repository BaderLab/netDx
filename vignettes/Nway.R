
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
idx <- union(which(pam50 %in% c("Normal-like","HER2-enriched")), 
	which(is.na(st2)))
idx <- union(idx, which(is.na(pam50)))
pID <- colData(brca)$patientID
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

pam50 <- colData(brca)$PAM50.mRNA
colData(brca)$pam_mod <- pam50

smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- brca[[1]][,notdup]


## ----eval=TRUE-----------------------------------------------------------
pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- gsub(" ",".",colData(brca)$pam_mod)

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

# Average accuracy
acc <- matrix(NA,ncol=length(st),nrow=numSplits) 
colnames(acc) <- st 

for (k in 1:numSplits) { 
	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
	oF <- sprintf("%s/Split%i_predictionResults.txt",outDir,k)
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",sprintf("%s_SCORE",st))]
	write.table(tmp,file=oF,sep="\t",col=TRUE,row=FALSE,quote=FALSE)

	# label-specific accuracy
	ctr <- 1
	for (m in st) { 
		tmp2 <- subset(tmp, STATUS==m)
		acc[k,ctr] <- (sum(tmp2$STATUS==tmp2$PRED_CLASS)/nrow(tmp2))*100
		ctr <- ctr+1
	}
}
print(round(acc,2))

print("Confusion matrix")
res <- out$Split1$predictions
print(table(res[,c("STATUS","PRED_CLASS")]))



## ------------------------------------------------------------------------
sessionInfo()

