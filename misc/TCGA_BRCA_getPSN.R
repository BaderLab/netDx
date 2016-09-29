#' generate integrated PSN for TCGA breast cancer data
rm(list=ls())

# patients
geneFile	<- "/home/spai/tmp/TCGA_BRCA_freeze160826/LumA/tmp/GENES.txt"
# net names
netFile		<- "/home/spai/tmp/TCGA_BRCA_freeze160826/LumA/tmp/NETWORKS.txt"
# binary nets
netDir		<- "/home/spai/tmp/TCGA_BRCA_freeze160826/LumA/tmp/INTERACTIONS"
# we are going to take union of FS pathways for each class so we need
# the pathway scores for each class
pTFile_LumA	<- "/home/spai/tmp/TCGA_BRCA_freeze160826/LumA/GM_results/LumA_pathway_CV_score.txt"
pTFile_other<- "/home/spai/tmp/TCGA_BRCA_freeze160826/other/GM_results/other_pathway_CV_score.txt"
# where the integrated net will be written
outDir		<- "/home/spai/tmp/TCGA_BRCA_freeze160826/LumA/output_nets_score10"

# -------
if (file.exists(outDir)) unlink(outDir,recursive=TRUE)
dir.create(outDir)

require(netDx)

cutoff <-10 
corrThresh <-0.7 
aggFun <- "MAX"

# run for LumA
cat("\nLumA\n")
pTally <- read.delim(pTFile_LumA,sep="\t",h=T,as.is=T)
pTally <- subset(pTally, pTally[,2]>=cutoff)
cat(sprintf("%i networks with score >=%i\n",nrow(pTally),cutoff))
colnames(pTally) <- c("NETWORK","WEIGHT")
aggNetFile <- netDx::writeWeightedNets(geneFile,netFile,netDir,
				pTally,outDir,
				filterEdgeWt=corrThresh,limitToTop=25,
				writeAggNet=aggFun,verbose=TRUE)

#aggNetFile <- sprintf("%s/aggregateNet_filterEdgeWt%1.2f_%s.txt",
#					  outDir,corrThresh,aggFun)
aggNet<- read.delim(aggNetFile,sep="\t",h=T,as.is=T)[,1:3]
colnames(aggNet)[3] <- "weight"

# compute shortest path within/across clusters
require(netDx.examples)
data(TCGA_BRCA)
pheno$STATUS[which(!pheno$STATUS %in% "LumA")] <- "other"
colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity
x <- compareShortestPath(aggNet, pheno)
write.table(pheno,file=sprintf("%s/LumA_pheno.txt",outDir),sep="\t",
			col=T,row=F,quote=F)

rm(pTally)

###warning("nets selected in both classes will be overwritten by second class' weight\n")
#### now run for other
###cat("\nOther\n")
###pTally <- read.delim(pTFile_other,sep="\t",h=T,as.is=T)
###pTally <- subset(pTally, pTally[,2]>=cutoff)
###cat(sprintf("%i networks with score >=%i\n",nrow(pTally),cutoff))
###colnames(pTally) <- c("NETWORK","WEIGHT")
###writeWeightedNets(geneFile,netFile,netDir,pTally,outDir,
###				  filterEdgeWt=corrThresh,
###				  writeMeanNet=FALSE,verbose=FALSE)
