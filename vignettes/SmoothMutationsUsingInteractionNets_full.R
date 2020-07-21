## ----eval=TRUE----------------------------------------------------------------



## ----eval=TRUE----------------------------------------------------------------
set.seed(8)
suppressWarnings(suppressMessages(require(netDx)))
suppressWarnings(suppressMessages(require(MultiAssayExperiment)))


## ----eval=TRUE----------------------------------------------------------------
genoFile <- paste(path.package("netDx"),"extdata",
                  "OV_mutSmooth_geno.txt",
	sep=getFileSep())
print(genoFile)
geno <- read.delim(genoFile,sep="\t",header=TRUE,as.is=TRUE)

phenoFile <- paste(path.package("netDx"),"extdata",
                   "OV_mutSmooth_pheno.txt",
                   sep=getFileSep())
pheno <- read.delim(phenoFile,sep="\t",header=TRUE,as.is=TRUE)


## ----eval=TRUE----------------------------------------------------------------
colnames(geno) <- gsub("-",".",colnames(geno))
pheno$ID <- gsub("-",".",pheno$ID)
rownames(pheno) <- pheno$ID


## ----eval=TRUE----------------------------------------------------------------
netFile <- paste(path.package("netDx"),"extdata",
                 "CancerNets.txt",
	sep=getFileSep())
cancerNets <- read.delim(netFile,sep="\t",header=T,as.is=T)
head(cancerNets)


## ----eval=TRUE----------------------------------------------------------------
message("* Running label prop")
require(doParallel)
# Start the node clusters for parallel propagation
cl <- makeCluster(1L)
registerDoParallel(cl)
prop_net <- smoothMutations_LabelProp(geno,cancerNets,cl,
                                      no_cores=10L)
stopCluster(cl)


## ----eval=TRUE----------------------------------------------------------------
#Set the name of the project and future resulting directory
#Apply binarization
genoP <- thresholdSmoothedMutations(
   prop_net,geno,"OV_CancerNets"
   )


## ----eval=TRUE----------------------------------------------------------------
#Setup to build the predictor
pathwayList <- readPathways(
   fetchPathwayDefinitions("January",2018)
   )
# limit num pathways to speed example
##pathwayList <- pathwayList[sample(1:length(pathwayList),100,replace=FALSE)]
exprdat <- SummarizedExperiment(genoP, colData=pheno)
objList <- list(genetic=exprdat)


## ---- eval=TRUE---------------------------------------------------------------
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


## ----eval=TRUE----------------------------------------------------------------
exprdat <- SummarizedExperiment(genoP, colData=pheno)
objList <- list(genetic=exprdat)

groupList <- list()
groupList[["genetic"]] <- pathwayList #names for groupList and objList now match

dataList <- MultiAssayExperiment(objList,pheno)


## ----eval=TRUE----------------------------------------------------------------
outDir <- tempdir()
#Run the predictor as usual
out <- buildPredictor(dataList=dataList,groupList=groupList,
  makeNetFunc=makeNets, ### custom network creation function
  outDir=outDir, ## absolute path
  numCores=10L, featScoreMax=10L, featSelCutoff=9L,
  numSplits=10L
)


## ----eval=FALSE---------------------------------------------------------------
## numSplits <- 2
## st <- unique(colData(dataList)$STATUS)
## acc <- c()         # accuracy
## predList <- list() # prediction tables
## 
## featScores <- list() # feature scores per class
## for (cur in unique(st)) featScores[[cur]] <- list()
## 
## for (k in 1:numSplits) {
## 	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
## 	# predictions table
## 	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
## 	                 sprintf("%s_SCORE",st))]
## 	predList[[k]] <- tmp
## 	# accuracy
## 	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
## 	# feature scores
## 	for (cur in unique(st)) {
## 	   tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
## 	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
## 	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
## 	}
## }


## ----eval=FALSE---------------------------------------------------------------
## predPerf <- plotPerf(predList, predClasses=st)


## ----eval=FALSE---------------------------------------------------------------
## featScores2 <- lapply(featScores, getNetConsensus)
## summary(featScores2)
## 
## featSelNet <- lapply(featScores2, function(x) {
##     callFeatSel(x, fsCutoff=1, fsPctPass=0)
## })
## print(head(featScores2[["deceased"]]))


## ----eval=TRUE----------------------------------------------------------------
sessionInfo()

