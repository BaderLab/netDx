# BRCA example with nested cv

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

subtypes<- c("LumA")
pheno$STATUS[which(!pheno$STATUS %in% subtypes)] <- "other"
subtypes <- c(subtypes,"other") # add residual

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
    path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)
head(pathwayList)


BRCA_makeNets <- function(dataList, groupList, netDir,...) {
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

rootDir <- "/home/shraddhapai/BaderLab/2017_BRCA/output/"
dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/BRCA_%s",rootDir,dt)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=BRCA_makeNets, ### custom network creation function
   outDir=megaDir,
   numCores=10L,nFoldCV=10L, CVcutoff=9L,numSplits=25L)
