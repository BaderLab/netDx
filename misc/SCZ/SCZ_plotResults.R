#' plot BRCA results
rm(list=ls())
require(netDx)
require(netDx.examples)
rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2018_CommonMind"
inDir <- sprintf("%s/input",rootDir)
geneFile <- sprintf("%s/anno/gencode.v27lift37.annotation.gtf.geneidsymbols.txt",rootDir)
pathFile <- sprintf("%s/anno/Human_AllPathways_November_01_2017_symbol.gmt",
	rootDir)

genes <- read.delim(geneFile,sep="\t",h=T,as.is=T)
gene_id <- genes[,4]
dpos <- regexpr("\\.",gene_id)
gene_id <- substr(gene_id,1,dpos-1)
gene_name <- genes[,5]

options(StringsAsFactors=FALSE)

# prepare input
pheno <- read.delim(sprintf("%s/CMC_MSSM-Penn-Pitt_Clinical.csv",
	inDir),sep=",",h=T,as.is=T)
xpr <- read.delim(sprintf("%s/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_gene-adjustedSVA-dataNormalization-includeAncestry-adjustedLogCPM.tsv.gz",inDir),sep="\t",h=T,as.is=T)

not_in <- which(!xpr[,1] %in% gene_id)
cat(sprintf("%i ids not in dictionary; excluding\n", length(not_in)))
if (length(not_in)>0) xpr <- xpr[-not_in,]

midx <- match(xpr[,1],gene_id)
if (all.equal(gene_id[midx],xpr[,1])!=TRUE) {
	cat("don't match\n")
	browser()
}
xpr_genes <- gene_name[midx]

idx <- which(pheno$Dx %in% c("Control","SCZ"))
pheno <- pheno[idx,]
pheno$STATUS <- pheno$Dx
pheno$ID <- pheno$DLPFC_RNA_Sequencing_Sample_ID

# filter genes by those profiled
pathwayList <- readPathways(pathFile)
for (k in names(pathwayList)) 
	pathwayList[[k]] <- intersect(pathwayList[[k]],xpr_genes)

inDir <- sprintf("%s/output/SCZ_180131/pred",rootDir)
outDir <- sprintf("%s/output/SCZ_180131/plot",rootDir)
if (!file.exists(outDir)) dir.create(outDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf <- plotPerf(inDir, predClasses=predClasses)
dev.off()

auc <- unlist(lapply(predPerf,function(x) x$auroc))
aupr <- unlist(lapply(predPerf,function(x) x$aupr))
acc <- unlist(lapply(predPerf,function(x) x$accuracy))
cat(sprintf("Perf: AUROC=%1.2f ; AUPR = %1.2f ; acc = %1.2f%%\n",
	mean(auc),mean(aupr),mean(acc)))

featScores <- getFeatureScores(inDir,predClasses=predClasses)
featSelNet <- lapply(featScores, function(x) {
    callFeatSel(x, fsCutoff=10, fsPctPass=0.9)
})


netInfoFile <- sprintf("%s/inputNets.txt",inDir)
netInfo <- read.delim(netInfoFile,sep="\t",h=FALSE,as.is=TRUE)
EMap_input <- writeEMapInput_many(featScores,pathwayList,
      netInfo,outDir=outDir)
