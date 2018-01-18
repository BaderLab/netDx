# Ependymoma
rm(list=ls())

require(netDx)
require(netDx.examples)

rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"
inDir <- sprintf("%s/input",rootDir)
outDir <- sprintf("%s/output",rootDir)

xpr <- read.delim(sprintf("%s/original_data/Toronto-comparison2-without-spinals/TOR-ST-PFPURE-PFMIX-SEP16.gct",inDir),skip=2,h=T,as.is=T)
rownames(xpr) <- xpr[,1]
xpr <- xpr[,-(1:2)]
sampType <- scan(sprintf("%s/original_data/Toronto-comparison2-without-spinals/TOR-ST-PFPURE-PFMIX-SEP16.cls",inDir),skip=2)
sampType <- as.integer(sampType)

# from Ruth
#  st = 0, PFPURE = 1 and PFMIX = 2
pheno <- data.frame(ID=colnames(xpr),INT_STATUS=sampType)
pheno$ID <- as.character(pheno$ID)
st <- c("ST","PFPURE","PFMIX")
pheno$STATUS <- st[sampType+1]

# exclude ST
idx <- which(pheno$STATUS=="ST") 
pheno <- pheno[-idx,]
xpr <- xpr[,-idx]
xpr <- log(xpr+1)

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
    path.package("netDx.examples"))
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
megaDir <- sprintf("%s/Epen_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=10L,nFoldCV=3L, CVcutoff=2L,numSplits=10L)
