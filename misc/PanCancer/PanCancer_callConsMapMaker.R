#' generate consensus maps for all PanCancer feature selection
#' predictors
require(netDx)
require(netDx.examples)

# root directory where all pancancer information is stored on the running
# machine
dirBase <- "/Users/shraddhapai/Documents/Research/BaderLab"

# with feature selection
dirList <- list(
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/output/featSel_incMutRPPA_170223",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222",dirBase)
	)

# directories with input files
inDirs <- list( 
	GBM=sprintf("%s/2017_TCGA_GBM/input",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/input",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/input",dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/input",dirBase)
)

# common data sources
# RNA-based pathway nets
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# ------------------------------------------------------------------------
### KIRC
### load genes in this dataset, to create an input set that only has genes
### measured here.

# somatic mutations
mutFile		<- sprintf("%s/KIRC_core_somatic_mutations.txt",inDirs$KIRC)
mut_genes <- system(sprintf("cut -f1 %s", mutFile),intern=TRUE)[-1]
mutList <- pathwayList; 
names(mutList) <- paste("MUT_",names(mutList),sep="")
mutList <- lapply(mutList, function(x) x[which(x %in% mut_genes)])
ln <- unlist(lapply(mutList, length))
if (any(ln<1)) mutList <- mutList[-which(ln<1)]

# gene expression
xprFile <- sprintf("%s/KIRC_mRNA_core.txt",inDirs$KIRC)
xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
bpos <- regexpr("\\|", xpr_genes)
xpr_genes <- substr(xpr_genes, 1,bpos-1)
xprList <- pathwayList
xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])

netInfo <- list(
	RNA=list(pattern=".profile$", netSet=xprList),
	mutations=list(pattern="MUT_", netSet=mutList),
	clinical=list(pattern="_cont",netSet="base")
	)

cat("got here\n")

source("PanCancer_featSel_writeConsensusMap.R")
writeConsensusMap(datDir=dirList$KIRC, predClasses=c("SURVIVEYES","SURVIVENO"),
	netInfo=netInfo)
									


rm(netInfo)
