# get nets sampled in design d to see if themes from real data are identified

require(netDx)
require(netDx.examples)

rootDir <- "/home/shraddhapai/BaderLab"
indir <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathRandom_designD_170713"

pathFile <- sprintf("%s/extdata/Human_160124_AllPathways.gmt",
   path.package("netDx.examples"))
pathwayList <- readPathways(pathFile)

nets <- list()
nets[["YES"]] <- list()
nets[["NO"]] <- list()
for (k in 1:100) {
	system(sprintf("cd %s/rng%i/SURVIVEYES/tmp; ls *.profile > %s/nets/rng%i_YES_nets.txt",indir,k,indir,k))
	system(sprintf("cd %s/rng%i/SURVIVENO/tmp; ls *.profile > %s/nets/rng%i_NO_nets.txt",indir,k,indir,k))
	
	dat <- read.delim(sprintf("%s/nets/rng%i_YES_nets.txt",
		indir,k),sep="\t",h=F,as.is=T)
	nets[["YES"]][[k]] <- dat[,1]
	dat <- read.delim(sprintf("%s/nets/rng%i_NO_nets.txt",
		indir,k),sep="\t",h=F,as.is=T)
	nets[["NO"]][[k]] <- dat[,1]
}

rnaFile <- sprintf("%s/PanCancer_KIRC/input/KIRC_mRNA_core.txt",
	rootDir)
rna <- read.delim(rnaFile,sep="\t",h=T,as.is=T)
genes <- sub("mRNA_","",colnames(rna))
dpos <- regexpr("\\.",genes)
genes <- substr(genes,1,dpos-1)[-1]
rm(rna)

# make matrix
for (g in c("YES","NO")) {
	cat(sprintf("Group %s\n",g))
	allp <- unique(unlist(nets[[g]]))
	mat <- matrix(0,nrow=length(allp),ncol=100)
	rownames(mat) <- allp 
	for (k in 1:100) {
		idx <- which(rownames(mat) %in% nets[[g]][[k]])
		mat[idx,k] <- 1
	}

	# write gmt
	cat("\tgmt\n")
	outFile <- sprintf("%s/nets/%s.gmt",indir,g)
	system(sprintf("cat /dev/null > %s",outFile))
	netScores <- rowSums(mat)
print(summary(netScores))
	rownames(mat) <- sub(".profile","",rownames(mat))
	mat <- mat[which(netScores>=7),]
	netScores <- netScores[which(netScores>=7)]
	for (curp in rownames(mat)) {
		curgenes <- intersect(pathwayList[[curp]],genes)
		cat(sprintf("%s\t%s\t%s\n",curp,curp,
			paste(curgenes,collapse="\t")),
			file=outFile,append=TRUE) 
	}
	
	# write tally
	outFile <- sprintf("%s/nets/%s_attr.gmt",indir,g)
	df <- data.frame(netName=rownames(mat),score=netScores)
	write.table(df,file=outFile,sep="\t",col=T,row=F,quote=F)
}

