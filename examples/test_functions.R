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
