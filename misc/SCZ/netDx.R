# Ependymoma
rm(list=ls())

require(netDx)
require(netDx.examples)

rootDir <- "/home/shraddhapai/BaderLab/2018_CommonMind"
inDir <- sprintf("%s/input",rootDir)
outDir <- sprintf("%s/output",rootDir)
pathFile <-sprintf("%s/anno/Human_AllPathways_November_01_2017_symbol.gmt",
	rootDir)

geneFile <- sprintf("%s/anno/gencode.v27lift37.annotation.gtf.geneids_chroms.txt.gz",rootDir)
genes <- read.delim(geneFile,sep="\t",h=F,as.is=T)
nm <- strsplit(genes[,5]," ")
gene_id <- unlist(lapply(nm,function(x)x[[1]]))
dpos <- regexpr("\\.",gene_id)
gene_id <- substr(gene_id,1,dpos-1)
gene_name <- unlist(lapply(nm,function(x) x[[2]]))

options(StringsAsFactors=FALSE)
pheno <- read.delim(sprintf("%s/CMC_MSSM-Penn-Pitt_Clinical.csv",
	inDir),sep=",",h=T,as.is=T)
xpr <- read.delim(sprintf("%s/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv.gz",inDir),sep="\t",h=T,as.is=T,nrow=1000)

not_in <- which(!xpr[,1] %in% gene_id)
cat(sprintf("%i ids not in dictionary; excluding\n", length(not_in)))
if (length(not_in)>0) xpr <- xpr[-not_in,]

midx <- match(xpr[,1],gene_id)
if (all.equal(gene_id[midx],xpr[,1])!=TRUE) {
	cat("don't match\n")
	browser()
}
xpr_genes <- gene_name[midx]
rownames(xpr)<- xpr_genes; rm(gene_id,gene_name)
xpr <- xpr[,-1]
common <- intersect(colnames(xpr),pheno$DLPFC_RNA_Sequencing_Sample_ID)
xpr <- xpr[,which(colnames(xpr) %in% common)]
pheno <- pheno[which(pheno$DLPFC_RNA_Sequencing_Sample_ID %in% common),]

midx <- match(pheno$DLPFC_RNA_Sequencing_Sample_ID,colnames(xpr))
if (all.equal(colnames(xpr)[midx],pheno$DLPFC_RNA_Sequencing_Sample_ID)!=TRUE) {
	cat("pheno don't match"); browser()
}
xpr <- xpr[,midx]

# limit by dx
idx <- which(pheno$Dx %in% c("Control","SCZ"))
xpr <- xpr[,idx]
pheno <- pheno[idx,]
pheno$STATUS <- pheno$Dx
pheno$ID <- pheno$DLPFC_RNA_Sequencing_Sample_ID
    
pathwayList <- readPathways(pathFile)
head(pathwayList)

makeNets <- function(dataList, groupList, netDir,...) {
	netList <- c()
	# make RNA nets: group by pathway
	if (!is.null(groupList[["rna"]])) { 
	netList <- makePSN_NamedMatrix(dataList$rna, 
					rownames(dataList$rna),
			   	groupList[["rna"]],netDir,verbose=FALSE, 
			  	writeProfiles=TRUE,...) 
	netList <- unlist(netList)
	cat(sprintf("Made %i RNA pathway nets\n", length(netList)))
	}
	cat(sprintf("Total of %i nets\n", length(netList)))
	return(netList)
}

dt <- format(Sys.Date(),"%y%m%d")
megaDir <- sprintf("%s/SCZ_%s",outDir,dt)
if (!file.exists(megaDir)) dir.create(megaDir)

gps <- list(rna=pathwayList)
dats <- list(rna=xpr)

browser()

runPredictor_nestedCV(pheno,
   dataList=dats,groupList=gps,
   makeNetFunc=makeNets, ### custom network creation function
   outDir=sprintf("%s/pred",megaDir),
   numCores=4L,nFoldCV=10L, CVcutoff=9L,numSplits=10L)
