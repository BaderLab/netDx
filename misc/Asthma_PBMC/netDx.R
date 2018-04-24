# Ependymoma
rm(list=ls())

require(GEOquery)
require(netDx)
require(netDx.examples)
require(org.Hs.eg.db)

rootDir <- "/home/shraddhapai/BaderLab/2018_AsthmaPBMC"
inDir <- sprintf("%s/input",rootDir)
outDir <- sprintf("%s/output",rootDir)
pathFile <-sprintf("%s/anno/Human_AllPathways_February_01_2018_symbol.gmt",
	rootDir)

dat <- getGEO(filename=sprintf("%s/GSE40732_series_matrix.txt.gz",inDir),
	GSEMatrix=TRUE)
xpr <- exprs(dat)
pheno <- pData(dat)
# map GB ID to symbol
x <- mapIds(org.Hs.eg.db, keys=rownames(xpr), column="SYMBOL",keytype="ACCNUM",
	multiVals="first")
common <- intersect(names(x),rownames(xpr))
xpr <- xpr[which(rownames(xpr) %in% common),]
x <- x[which(names(x) %in% common)]

midx <- match(rownames(xpr),names(x))
gnames <- x[midx]
agg <- aggregate(xpr, by=list(gene_name=gnames),FUN=mean)
xpr <- agg[,-1]
rownames(xpr) <- agg[,1]

pheno <- pheno[,c("geo_accession","characteristics_ch1")]
st <- rep(NA, nrow(pheno))
st[which(pheno[,2] %in% "asthma: FALSE")] <- "control"
st[which(pheno[,2] %in% "asthma: TRUE")] <- "asthma"
pheno[,2] <- st
colnames(pheno) <- c("ID","STATUS")
pheno[,1] <- as.character(pheno[,1])

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
megaDir <- sprintf("%s/basic_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=4L,nFoldCV=10L, CVcutoff=9L,numSplits=100L,CVmemory=13L)	
