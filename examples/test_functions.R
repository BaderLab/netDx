#' Example snippets that should work in contained manner.
rm(list=ls())

if ("package:netDx" %in% search()) detach(package:netDx,unload=TRUE);
require(netDx)

tmpDir <- "~/tmp/netDx_tests"

## makePSN_Named_Matrix
data(TCGA_mini,pathwayList); 
out <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList, tmpDir,writeProfiles=TRUE)
rm(xpr,pheno,cnv_GR,pathwayList)

## mapNamedRangesToSets
data(genes,pathwayList); 
gene_GR<- GRanges(genes$chrom,IRanges(genes$txStart,genes$txEnd),name=genes$name2)
path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
rm(path_GRList,pathwayList,gene_GR,genes)

## makePSN_RangeSets
data(pathway_GR,TCGA_mini); 
netList <- makePSN_RangeSets(cnv_GR,pathway_GR,tmpDir)
rm(pathway_GR,xpr,pheno,cnv_GR,netList)

## run_GeneMANIA
GM_db <- sprintf("%s/extdata/GM_db", path.package("netDx"))
GM_query <- sprintf("%s/extdata/GM_query.txt",path.package("netDx"))
runGeneMANIA(GM_db, GM_query,tmpDir)
rm(GM_db,GM_query)

## GM_parseReport
x <- runGeneMANIA(GM_db,GM_query,tmpDir)
GM_parseReport(x)

## readPathways
pathFile <- sprintf("%s/extdata/pathways.gmt",path.package("netDx"))
x <- readPathways(pathFile)
rm(x,pathFile)

## GM_createDB
data(TCGA_mini,pathwayList); 
n <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList,sprintf("%s/nets",tmpDir),
	writeProfiles=TRUE); 
db <- GM_createDB(sprintf("%s/nets",tmpDir),pheno$ID,tmpDir)
rm(xpr,pheno,cnv_GR,n,db)

## GM_getQueryROC
data(TCGA_mini); 
prankFile <- sprintf("%s/extdata/GM_PRANK.txt", 
					 path.package("netDx"))
x <- GM_getQueryROC(prankFile, pheno, "LumA")

## GM_networkTally
netDir <- sprintf("%s/extdata/GM_NRANK",path.package("netDx"))
netFiles <- sprintf("%s/%s", netDir,dir(netDir,pattern="NRANK$"))
pTally <- GM_networkTally(netFiles)

## GM_writeQueryFile
data(TCGA_mini); GM_writeQueryFile(pheno$ID[1:5], "all",nrow(pheno),
								   sprintf("%s/myquery.txt",tmpDir))
rm(xpr,pheno,cnv_GR)

## GM_OneVAll_getClass
data(predRes); predClass <- GM_OneVAll_getClass(predRes)
rm(predRes,predClass)

## makeCVqueries
data(TCGA_mini); x <- makeCVqueries(pheno$ID)
rm(xpr,pheno,cnv_GR,x)

## makeRandomNetworks
data(TCGA_mini); x <- makeRandomNetworks(pheno$ID, numNets=10L,
	outDir=tmpDir)
rm(xpr,pheno,cnv_GR,x)

## getSimilarity
data(TCGA_mini); x <- getSimilarity(xpr)
mySim <- function(x) cor(x,method="kendall")
x <- getSimilarity(xpr,customFunc=mySim)
rm(xpr,pheno,cnv_GR,x)

## splitTestTrain
data(TCGA_mini); x <- splitTestTrain(pheno,predClass="LumA")

## splitTestTrain_partition
data(TCGA_mini); x <- splitTestTrain_partition(pheno,predClass="LumA")

## perfCalc
data(confmat); x <- perfCalc(confmat)

## cleanPathwayName
cleanPathwayName("7-(3-AMINO-3-CARBOXYPROPYL)-WYOSINE BIOSYNTHESIS%HUMANCYC%PWY-7286")

## sparsifyNet
require(reshape2)
data(TCGA_mini)
x <- melt(cor(xpr)) 
y <- sparsifyNet(x,outFile="tmp.txt")
rm(xpr,pheno,cnv_GR,x,y)

## getEnr
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
x <- getEnr(netDir, npheno, "case",netGrep=".txt$")
print(x$orig_rat) 
rm(npheno,x,netDir)

## cliqueFilterNets
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
x <- cliqueFilterNets(netDir,npheno,"~/tmp",predClass="case",netGrep="txt$",numReps=500)
print(x)
rm(npheno,x,netDir)

## countIntType_batch
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
countIntType_batch(sprintf("%s/BOTH_EQUAL.txt", netDir),npheno[1:100,1],npheno[101:200,1])

## countPatientsInNet
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
x <- countPatientsInNet(netDir,dir(netDir,pattern="txt$"), npheno[,1])

## getOR
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
x <- countPatientsInNet(netDir,dir(netDir,pattern="txt$"), npheno[,1])
y <- getOR(x,npheno,"case",colnames(x)[1]) # should give large OR
y <- getOR(x,npheno,"case",colnames(x)[2]) # should give OR of 0
y <- getOR(x,npheno,"case",colnames(x)[3]) # should give OL of 1

## pruneNets
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
pruneNets(netDir,"~/tmp",filterIDs=npheno[1:10,],netSfx="txt$")

## updateNets
data(npheno)
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
netmat <- countPatientsInNet(netDir,dir(netDir,pattern="txt$"), npheno[,1])
updateNets(netmat, npheno,writeNewNets=FALSE)

## writeNetsSif
netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
netFiles <- sprintf("%s/%s", netDir, dir(netDir,pattern="txt$"))
writeNetsSIF(netFiles,"~/tmp/merged.sif",netSfx=".txt")

## CV_perf
data(TCGA_full)
pDir <- sprintf("%s/extdata/GM_PRANK",path.package("netDx"))
pranks <- sprintf("%s/%s",pDir,dir(pDir,pattern="PRANK$")) # example PRANK
x <- CV_perf(pranks,pheno_full, "LumA","~/tmp")
x <- CV_perf(pranks,pheno_full, "LumA",justGetError=TRUE)

## doSubsetSelection
### warning - takes a long time to run. Recommend running on a multi-core
### machine with numCores set accordingly
GM_db <- sprintf("%s/extdata/GM_db",path.package("netDx"))
data(MB_pheno)
nets <- c("DKK2_cont","EMX2_cont","TNC_cont","WIF1_cont")
doSubsetSelection(methodName="greedy.backward",incNets=nets, 
	queryPool=MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
	GM_db=GM_db,pheno_DF=MB.pheno, outDir="~/tmp",num2return=103L)

## GM_runCV_nested
data(MB_pheno)
GM_db <- sprintf("%s/extdata/GM_db",path.package("netDx"))
GM_runCV_nested(outDir="~/tmp",pheno_DF=MB.pheno, predictClass="WNT",
	trainID_pred=MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
	numTrainSamps=103L,GM_db=GM_db)

## GM_runCV_featureSet
data(MB_pheno)
GM_db <- sprintf("%s/extdata/GM_db",path.package("netDx"))
GM_runCV_featureSet(MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
	"~/tmp",GM_db,103L)

## fetch_NetUnits
data(TCGA_mini,pathway_GR,pathwayList)
x <- getRegionOL(cnv_GR,pathway_GR)
y <- fetch_NetUnits(x,pathwayList, names(pathwayList))
y <- fetch_NetUnits(x,pathwayList, names(pathwayList),trackMapping_detail=TRUE)

