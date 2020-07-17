
#Set working directory----
set.seed(8)

numCores <- 8L

#Load libraries----
#library("data.table")
#library("matrixStats")
#library("plyr")
#library("doParallel")
#library("parallel")
#library("scater")
library("netDx")
require(MultiAssayExperiment)

#Set input and output paths ----
#Input
#nets_path="input/nets_l.rda"
#datasets_path="input/datasets.rda"
#Output
outDir <- "/output"

#Load data
#load(nets_path)
#load(datasets_path)


# describe dataset here.
# These are data from XXXX, XXX samples [cite paper].

# describe data type - binary matrix, with rows corresponding to genes
genoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_geno.txt",
	sep=getFileSep())
geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)

phenoFile <- paste(path.package("netDx"),"extdata","OV_mutSmooth_pheno.txt",
	sep=getFileSep())
pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)

# R behaves strangely with hyphens in rownames. Pre-empt by replacing "-" with 
# "." for all sample names
colnames(geno) <- gsub("-",".",colnames(geno))
pheno$ID <- gsub("-",".",pheno$ID)
rownames(pheno) <- pheno$ID

netFile <- paste(path.package("netDx"),"extdata","CancerNets.txt",
	sep=getFileSep())
cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)

# remove genes not in network, from geno
message("* Excluding genes not present in interaction nets")
noData <- setdiff(rownames(geno),rownames(cancerNets))
if (length(noData)>0) {
	message(paste(length(noData), 
		" genes not present in cancer nets; excluding",sep=""))
	geno <- geno[-which(rownames(geno) %in% noData),]
}
#geno <- complete_m(geno,cancerNets)


message("* Running label prop")
require(doParallel)
cl <- makeCluster(numCores)
registerDoParallel(cl)
prop_net <- prop_m(geno,cancerNets,cl,no_cores=numCores)
stopCluster(cl)

#Set the name of the project and future resulting directory
#Apply binarization
genoP <- discr_prop_l(prop_net,geno,"OV_CancerNets")

#Setup to build the predictor
pathwayList <- readPathways(fetchPathwayDefinitions("January",2018))
exprdat <- SummarizedExperiment(genoP, colData=pheno)
objList <- list(genetic=exprdat)

dataList <- MultiAssayExperiment::MultiAssayExperiment(objList,pheno)

groupList <- list()
groupList[["genetic"]] <- pathwayList #names for groupList and objList now match

#PSN function passed to netDx ----
makeNets <- function(dataList,groupList,netDir,numCores,...) {
  netList <- c(); netList2 <- c()
  
  # create genetic nets
  if (!is.null(groupList[["genetic"]])) {
	netList <- makeMutNets(dataList[["genetic"]],
		groupList[["genetic"]],
		netDir,numC=numCores)
  }
  cat(sprintf("\t%i genetic-pathway nets\n", length(netList)))
  cat(sprintf("Total of %i nets\n", length(netList)))
  
  return(netList)
}

# g geno matrix, genes by patients (columns) - binary
# pList list of genesets
# outDir - dir where nets are to be written
makeMutNets <- function(g,pList,oDir,numC) {
  g <- t(g) # transpose to have genes as columns
  cl	<- makeCluster(numC)
  registerDoParallel(cl)
  
  numPat <- c()
  netList <- foreach(k=1:length(pList)) %do% {
    idx <- which(colnames(g) %in% pList[[k]])
    
    if (length(idx)>0) {
      has_mut <- rowSums(g[,idx,drop=FALSE])
      has_mutp <- names(has_mut)[which(has_mut>0)]
      
      if (length(has_mutp)>=6) {
        cat(sprintf("%s: %i patients\n", names(pList)[k],
                    length(has_mutp)))
        #numPat <- c(numPat, length(has_mutp))
        pat_pairs <- t(combinat::combn(has_mutp,2));
        pat_pairs <- cbind(pat_pairs,1);
        outFile <- sprintf("%s/%s_cont.txt",oDir,names(pList)[k])
        write.table(pat_pairs, file=outFile,sep="\t",
                    col=FALSE,row=FALSE,quote=FALSE)
        basename(outFile)
      } else NULL
    } else {
      NULL
    }
  }
  stopCluster(cl)
  unlist(netList)
}

#Run the predictor as usual
buildPredictor(dataList=dataList,groupList=groupList,
  makeNetFunc=makeNets, ### custom network creation function
  outDir=outDir, ## absolute path
  numCores=numCores, numSplits=2L,debugMode=TRUE
)

# plot results.
# Show the ROC curve as in the breast cancer example
# Show top-scoring features
