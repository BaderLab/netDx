# Ependymoma
rm(list=ls())

require(netDx)
require(netDx.examples)

rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"
inDir <- sprintf("%s/input/netDx_prepared",rootDir)
outDir <- sprintf("%s/output",rootDir)
pathFile <-sprintf("%s/anno/Human_AllPathways_February_01_2018_symbol.gmt",
	rootDir)
load(sprintf("%s/Ependymoma_cohortMerged_180125.Rdata",inDir))

# exclude ST
idx <- which(pheno$STATUS=="ST") 
pheno <- pheno[-idx,]
xpr <- xpr[,-idx]
    
pathwayList <- readPathways(pathFile)
head(pathwayList)

makeNets <- function(dataList, groupList, netDir,...) {
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
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/Epen_prunedPathway_0.001_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)
pheno$STATUS <- droplevels(pheno$STATUS)

#### -----------------------------------------------------
### BEGIN PRUNING CODE
# apply pruning to proteomic data 
curwd <- getwd()
setwd("../PanCancer")
source("LMprune.R")
source("runLM.R")
source("silh.R")
require(cluster)
setwd(curwd)
for (nm in setdiff(names(dats),"clinical")) {
print(nm)
	#if (nrow(dats[[nm]])>10000) topVar <- 50 else topVar <- 100
	pdf(sprintf("%s/%s_prune.pdf",megaDir,nm))
	prune <- LMprune(dats[[nm]],pheno$STATUS,topVar=100)
	dev.off()
	if (!is.na(prune)) {
		if (prune$bestThresh < 1) {
		res <- prune$res
		res <- subset(res, adj.P.Val < 0.001)
		tmp <- dats[[nm]];orig_ct <- nrow(tmp)
		tmp <- tmp[which(rownames(tmp)%in% rownames(res)),]
		dats[[nm]] <- tmp
		cat(sprintf("%s: Pruning with cutoff %1.2f\n", nm,prune$bestThresh))
		cat(sprintf("\t%i of %i left\n", nrow(tmp),orig_ct))
		}
	} else {
		cat(sprintf("%s: not pruning\n",nm))
	}
}
#### ----------------------------------------------------------

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=8L,nFoldCV=10L, CVcutoff=9L,numSplits=10L,startAt=1L)
