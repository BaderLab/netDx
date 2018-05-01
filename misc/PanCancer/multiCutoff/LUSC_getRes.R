#' plot GBM results with multiple CV cutoffs
rm(list=ls())
require(netDx)
require(reshape2)

mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output"
dirSet <- list(
	base="noPrune_180423",
	lasso="lasso_180426",
	lassoGenes="lassoGenes_180426",
	pamrGenes_sp2="pamrGenes_180427",
	pamrGenes_sp1="pamrGenes_sp1_180427"
)
settypes <- c("clinical","mir","rna","prot","cnv",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAcnv","all")

mega_auc <- list()
for (curdir in names(dirSet)) {
dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
	if (curdir %in% "base") rngMax <- 15
	else rngMax <- 20

auc_set <- list()
for (settype in settypes) {
	if (curdir %in% "lassoGenes") {
	rngDir <- paste(sprintf("%s/rng",dataDir), 3:rngMax,sep="")
	} else {
	rngDir <- paste(sprintf("%s/rng",dataDir), 1:rngMax,sep="")
	}


for (cutoff in 9) {
	c7 <- sprintf("%s/%s/cutoff%i/predictionResults.txt",
				  rngDir,settype,cutoff)
	torm <- c()
	for (idx in 1:length(c7)) {
		dat <- read.delim(c7[idx],sep="\t",h=T,as.is=T)
		x1 <- sum(dat$STATUS=="SURVIVEYES")
		x2 <- sum(dat$STATUS=="SURVIVENO")
		if (x1<1 & x2<1) torm <- c(torm, idx)
	}
	cat(sprintf("%i: removing %i\n", cutoff,length(torm)))
	if (length(torm)>0) c7 <- c7[-torm]
	postscript(sprintf("tmp.eps"))
	x <- plotPerf(c7,c("SURVIVEYES","SURVIVENO"))
	dev.off()
	y1 <- unlist(lapply(x,function(i) i$auroc))
	auc_set[[settype]] <- y1

	tmp <- c()
	cur <- auc_set[[settype]]
}
}
mega_auc[[curdir]] <- unlist(lapply(auc_set,mean))
}
dt <- format(Sys.Date(),"%y%m%d")
pdf(sprintf("LUSC_%s.pdf",dt),width=13,height=6); boxplot(mega_auc,main="LUSC",cex.axis=1.7,cex.main=2,las=1); dev.off()
