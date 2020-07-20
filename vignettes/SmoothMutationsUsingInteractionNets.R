# Set global enviroment variables ----
set.seed(8)
numCores <- 10L
numSplits <- 30L

# Load libraries ----
library("netDx")
require("MultiAssayExperiment")

# Set and load input data ----
# Ovarian serous cystadenocarcinoma (OV) data from Broad Institute TCGA Genome 
# Data Analysis Center (2016): 
# Firehose 2016_01_28 run. Broad Institute of MIT and Harvard. 
# doi:10.7908/C11G0KM9
# Data are rapresented in a binary matrix of genes as rows and 
# patients as columns
# Patients are divided in two classes based on their last follow up declared 
# status (172 deceased patients; 191 living)
# A gene is set to 1 for a patient if has a somatic mutation (
# single nucleotide polymorphism as well as di/tri/oligo-nucleotide 
# polymorphism) or indel 
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

# remove genes from geno matrix that are not in network
message("* Excluding genes not present in interaction nets")
noData <- setdiff(rownames(geno),rownames(cancerNets))
if (length(noData)>0) {
	message(paste(length(noData), 
		" genes not present in cancer nets; excluding",sep=""))
	geno <- geno[-which(rownames(geno) %in% noData),]
}

message("* Running label prop")
require(doParallel)
# Start the node clusters for parallel propagation
cl <- makeCluster(numCores)
registerDoParallel(cl)
prop_net <- smoothMutations_LabelProp(geno,cancerNets,cl,no_cores=numCores)
stopCluster(cl)

#Set the name of the project and future resulting directory
genoP <- thresholdSmoothedMutations(prop_net,geno,"OV_CancerNets")

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
out=buildPredictor(dataList=dataList,groupList=groupList,
  makeNetFunc=makeNets, ### custom network creation function
  outDir=outDir, ## absolute path
  numCores=numCores, numSplits=numSplits, debugMode=TRUE
)

#Derive the performances
st <- unique(colData(dataList)$STATUS)
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

#Plot the performances
predPerf <- plotPerf(predList, predClasses=st)

#Get top features and their score for deceased patient class
featScores2 <- lapply(featScores, getNetConsensus)
summary(featScores2)

featSelNet <- lapply(featScores2, function(x) {
    callFeatSel(x, fsCutoff=1, fsPctPass=0)
})
print(head(featScores2[["deceased"]]))

