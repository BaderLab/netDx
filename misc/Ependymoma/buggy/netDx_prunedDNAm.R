# Ependymoma - DNA methylation after Nature brain tumour paper
rm(list=ls())

require(netDx)
require(netDx.examples)

rootDir <- "/home/shraddhapai/BaderLab/2017_Ependymoma"
outDir <- sprintf("%s/output",rootDir)

phenoFile <- "/home/shraddhapai/BaderLab/2018_Epen_DNAm/input/GSE90496_pData.txt"
dnaFile <- "/home/shraddhapai/BaderLab/2018_Epen_DNAm/input/GSE90496_EPN_PFAB_beta.txt.gz"

# ----------------------
# input processing
pheno <- read.delim(phenoFile,sep="\t",h=T,as.is=T)
ttype <- pheno$characteristics_ch1
idx <- which(ttype %in% c("methylation class: EPN, PF A","methylation class: EPN, PF B"))
cat(sprintf("Got %i samples\n",length(idx)))
pheno <- pheno[idx,] # limit to EPN samples
cpos <- regexpr("sample", pheno$title)
bpos <- regexpr("\\[reference", pheno$title)
str <- as.integer(substr(pheno$title, cpos+7, bpos-2)) # get sample number
pheno$ID <- paste("SAMPLE", str,sep=".")
pheno <- pheno[,c("ID","characteristics_ch1")]
st <- rep("",nrow(pheno))
st[grep("PF A", pheno[,2])] <- "PFA"
st[grep("PF B", pheno[,2])] <- "PFB"
pheno$STATUS <- st

# ----------------------
makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["dnam"]])) { 
	netList <- makePSN_NamedMatrix(dataList$dnam, 
					rownames(dataList$dnam),
			   	groupList[["dnam"]],netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	netList <- unlist(netList)
	cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
	}
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/Epen_prunedOneNet_0.001_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

xpr <- read.delim(dnaFile,sep="\t",h=T,as.is=T)
rownames(xpr) <- paste("probe",1:nrow(xpr),sep="")
# match pheno and methylation values
midx <- match(colnames(xpr),pheno$ID)
if (all.equal(colnames(xpr),pheno$ID[midx])!=TRUE){
	cat("don't match\n")
browser()
}
pheno <- pheno[midx,]

gps <- list(dnam=list(dnam=rownames(xpr)))
dats <- list(dnam=xpr)

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
	if (nrow(dats[[nm]])>10000) topVar <- 50 else topVar <- 100
	pdf(sprintf("%s/%s_prune.pdf",megaDir,nm))
	prune <- LMprune(dats[[nm]],pheno$STATUS,topVar=topVar)
	dev.off()
	if (!is.na(prune)) {
		if (prune$bestThresh < 1) {
		res <- prune$res
		res <- subset(res, adj.P.Val < prune$bestThresh)
		
		require(dataExplore)
		cat("add hclust\n"); browser()
		
		tmp <- dats[[nm]];orig_ct <- nrow(tmp)
		tmp <- tmp[which(rownames(tmp)%in% rownames(res)),]
		dats[[nm]] <- tmp
		gps[[nm]] <- list(dnam=rownames(tmp))
		cat(sprintf("%s: Pruning with cutoff %1.2f\n", 
			nm,prune$bestThresh))
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
