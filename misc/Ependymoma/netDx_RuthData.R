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
megaDir <- sprintf("%s/Epen_2_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

pheno$STATUS <- droplevels(pheno$STATUS)

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=8L,nFoldCV=10L, CVcutoff=9L,numSplits=10L,startAt=8L)
