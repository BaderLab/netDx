#' call code that creates enrichment map.

setInfoFile <- "/Users/shraddhapai/Documents/Software/netDx/misc/PanCancer/featSel_pathways/KIRCpathway_locations.txt"
setInfo <- read.delim(setInfoFile,sep="\t",h=T,as.is=T)

#setInfo <- subset(setInfo, name %in% "pathOnly")
setInfo <- subset(setInfo, name %in% "pathOnly_noFS") #c("pathOnly80","pathOnly90","pathOnly95"))

# pathways only
datRoot <- "/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output"
outRoot <- "/Users/shraddhapai/DropBox/netDx/BaderLab/2017_PanCancer_Survival"

inDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/input"
xprFile <- sprintf("%s/KIRC_mRNA_core.txt",inDir)

# RNA-based pathway nets
require(netDx.examples)
require(netDx)
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
           path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# read genes measured in this set so we can keep ust those in EMap
xpr_genes <- scan(xprFile,nlines=1,what="character",quiet=TRUE)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes)
bpos <- regexpr("\\|", xpr_genes)
xpr_genes <- substr(xpr_genes, 1,bpos-1)
xprList <- pathwayList
xprList <- lapply(xprList, function(x) x[which(x %in% xpr_genes)])

source("writeConsensusNets_oneSet.R")
source("writeEMap.R")

for (curSet in 1:nrow(setInfo)) {
	cat(sprintf("%s\n",setInfo$name[curSet]))
	# write net scores
	od <- sprintf("%s/%s",outRoot,setInfo$outdir[curSet])
	if (!file.exists(od)) dir.create(od)

	netDir		<- sprintf("%s/netScores",od)
	dir.create(netDir)
	datDir		<- sprintf("%s/%s",datRoot, setInfo$dataDir[curSet])
	
	scorePfx	<- writeConsensusNets(datDir=datDir,
		outPfx=sprintf("%s/%s",netDir,setInfo$dataDir[curSet]),
		consCutoff=10,pctPass=0.7)

	# write gmt for enrichment map
	cat("Writing enrichment map\n")
	emapDir <- sprintf("%s/%s/EMap",outRoot,setInfo$outdir[curSet])
	dir.create(emapDir)
	for (gp in c("SURVIVEYES","SURVIVENO")) {
		cat(sprintf("\t%s\n",gp))
		nFile <- sprintf("%s_%s_netScores.txt",scorePfx, gp)
		writeEMap(nFile, xprList,netInfo=NULL,
				minScore=7,outPfx=sprintf("%s/%s",emapDir,gp))
	}
}
