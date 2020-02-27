suppressWarnings(suppressMessages(require(netDx)))
suppressMessages(library(curatedTCGAData))

set.seed(42) # make results reproducible

curatedTCGAData(diseaseCode="BRCA", assays="*",dry.run=TRUE)
brca <- suppressMessages(curatedTCGAData("BRCA",c("mRNAArray"),FALSE))

staget <- sub("[abcd]","",sub("t","",colData(brca)$pathology_T_stage))
staget <- suppressWarnings(as.integer(staget))
colData(brca)$STAGE <- staget

pam50 <- colData(brca)$PAM50.mRNA
pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
pam50[which(pam50 %in% "Luminal A")] <- "LumA"
colData(brca)$pam_mod <- pam50

tmp <- colData(brca)$PAM50.mRNA
idx <- union(which(tmp %in% c("Normal-like","Luminal B","HER2-enriched")),
		which(is.na(staget)))
pID <- colData(brca)$patientID
tokeep <- setdiff(pID, pID[idx])
brca <- brca[,tokeep,]

# remove duplicate asssays mapped to the same sample
smp <- sampleMap(brca)
samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
notdup <- samps[which(!duplicated(samps$primary)),"colname"]
brca[[1]] <- suppressMessages(brca[[1]][,notdup])

pID <- colData(brca)$patientID
colData(brca)$ID <- pID
colData(brca)$STATUS <- colData(brca)$pam_mod

groupList <- list()

# genes in mRNA data are grouped by pathways
pathList <- readPathways(
	fetchPathwayDefinitions("January",2018))
idx <- sample(1:length(pathList),10,F)
groupList[["BRCA_mRNAArray-20160128"]] <- pathList[idx]
# clinical data is not grouped; each variable is 
# its own feature
groupList[["clinical"]] <- list(
     age="patient.age_at_initial_pathologic_diagnosis",
	   stage="STAGE"
)

makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c() # initialize before is.null() check
	# (Pearson correlation)
	if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) { 
	netList <- makePSN_NamedMatrix(
				dataList[["BRCA_mRNAArray-20160128"]],
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
		simMetric="custom",customFunc=normDiff, 
		writeProfiles=FALSE,
		sparsify=TRUE,verbose=TRUE,...)
	}
	netList <- c(unlist(netList),unlist(netList2))
	return(netList)
}

outDir <- sprintf("%s/pred_output",tempdir()) # location for intermediate work
# set keepAllData=TRUE to not delete at the end of the predictor run.
# This can be useful for debugging.
numSplits <- 3L
out <- buildPredictor(
	dataList=brca,groupList=groupList,
  makeNetFunc=makeNets,outDir=outDir,
  numSplits=numSplits,featScoreMax=10L, featSelCutoff=9L,
	numCores=1L)

st <- unique(colData(brca)$STATUS)
acc <- c()         # accuracy
predList <- list() # prediction tables

featScores <- list() # feature scores per class
for (cur in unique(st)) featScores[[cur]] <- list()

for (k in 1:numSplits) { 
	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
	# predictions table
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
	                 sprintf("%s_SCORE",st))]
	predList[[k]] <- tmp 
	# accuracy
	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
	# feature scores
	for (cur in unique(st)) {
	   tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
	}
}

featScores2 <- lapply(featScores, getNetConsensus)
featSelNet <- lapply(featScores2, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=1)
})
print(featSelNet)
browser()

topPath <- gsub(".profile","",unique(unlist(featSelNet[["notLumA"]])))
topPath <- gsub("_cont.txt","",topPath)

# limit features to top ones
g2 <- list();
for (nm in names(groupList)) {
	cur <- groupList[[nm]]
	idx <- which(names(cur) %in% topPath)
	message(sprintf("%s: %i pathways", nm, length(idx)))
	if (length(idx)>0) g2[[nm]] <- cur[idx]
}

netDir <- sprintf("%s/final",outDir)
dir.create(netDir)
dir.create(sprintf("%s/profiles",netDir))

dat <- dataList2List(brca)
pheno <- dat$pheno
pheno_id <- setupFeatureDB(pheno,netDir)
createPSN_MultiData(dataList=dat$assays,groupList=g2,
			pheno=pheno_id,
			netDir=netDir,customFunc=makeNets,numCores=1,
			verbose=FALSE)
convertProfileToNetworks(
		netDir=sprintf("%s/profiles",netDir),
		outDir=sprintf("%s/INTERACTIONS",netDir),
)
#### rename
###networks <- read.delim(sprintf("%s/NETWORKS.txt",netDir),
###	sep="\t",header=FALSE,as.is=TRUE)
###networks <- networks[grep("profile$", networks[,2]),]
###networks[,2] <- sub(".profile","",networks[,2])
###
###for (i in 1:nrow(networks)) {
###	netid <- networks[i,1]
###	file.rename(from=sprintf("%s/INTERACTIONS/%s.txt",netDir,networks[i,2]),
###						to=sprintf("%s/INTERACTIONS/1.%i.txt",netDir,netid))
###}

###source("plotIntegratedPSN.R")
###source("writeWeightedNets.R")
###source("pruneNetByStrongest.R")
###suppressMessages(require(igraph))
###require(RColorBrewer)
###require(RCy3)
pheno_id <- pheno_id[,c("ID","STATUS")]
plotIntegratedPSN("LuminalA",pheno=pheno_id,netDir,topX=0.2)

