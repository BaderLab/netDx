# tally the frequency of strictNoFS sampled pathways, and make an emap.
require(netDx)
require(netDx.examples)

netList <- c()
for (f in dir(pattern="sampledNets")){
	dat <- read.delim(f,sep="\t",h=F,as.is=T)
	netList <- c(netList, basename(dat[,1]))
}

# count num times pathway shows up in noFS rounds
tally <- table(netList)
tally_int <- as.integer(tally)
names(tally_int) <- sub(".profile","",names(tally))
tally_int <- sort(tally_int,decreasing=TRUE)
tally <- data.frame(netName=names(tally_int), netScore=tally_int)
tally <- subset(tally, netScore>=10)

# load pathway and gene expression data
pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt", 
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

# limit to genes 
rnaFile <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/input/KIRC_mRNA_core.txt"
cat("Reading rna\n")
rna <- read.delim(rnaFile,sep="\t",h=T,as.is=T)
xpr_genes <- colnames(rna)[-1]
xpr_genes <- sub("mRNA_","",xpr_genes);
dpos <- regexpr("\\.",xpr_genes)
xpr_genes <- substr(xpr_genes,1,dpos-1)
rm(rna)
cat(sprintf("Set measured %i genes\n", length(xpr_genes)))

netAttrFile <- "strictNoFS_attr.txt"
outFile <- "strictNoFS.gmt"
if (file.exists(outFile)) unlink(outFile)
system(sprintf("touch %s",outFile))

# write node attributes
write.table(tally,file=netAttrFile,sep="\t",col=T,row=F,quote=F)
# write gmt
for (cur in tally$netName) {
	curGenes <- intersect(pathwayList[[cur]],xpr_genes)
	cat(sprintf("%s\t%s\t%s\n", cur,cur,
	paste(curGenes,collapse="\t")),file=outFile,append=TRUE)
}


