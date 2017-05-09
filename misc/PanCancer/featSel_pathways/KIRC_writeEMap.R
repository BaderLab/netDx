#' call code that creates enrichment map.

# pathways only
datDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/pathwaysOnly_170502"
inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/input"
scorePfx <- sprintf("%s/KIRC_thresh10_pctPass0.70",datDir)
netInfo <- sprintf("%s/inputNets.txt",datDir)
xprFile <- sprintf("%s/KIRC_mRNA_core.txt",inDir)

# RNA-based pathway nets
require(netDx.examples)
require(netDx)
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# gene expression
xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
bpos <- regexpr("\\|", xpr_genes)
xpr_genes <- substr(xpr_genes, 1,bpos-1)
xprList <- pathwayList
xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])

netInfo <- read.delim(netInfo,sep="\t",h=FALSE,as.is=T)

source("writeEMap.R")
for (gp in c("SURVIVEYES","SURVIVENO")) {
	nFile <- sprintf("%s_%s_netScores.txt",scorePfx, gp)
	writeEMap(nFile, xprList,netInfo=netInfo,
			outPfx=sprintf("%s_%s",datDir,gp))
	
}
