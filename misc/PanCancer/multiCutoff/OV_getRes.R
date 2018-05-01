#' plot GBM results with multiple CV cutoffs
rm(list=ls())
require(netDx)
require(reshape2)

mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/OV/output"

maxRng <- 100
settypes <- c("clinical","mir","rna","prot","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAdnam",
	"clinicalAcnv","all")
dirSet <- list(
	base="noPrune_180423",
	prune="pruneTrain_180419",
	lasso="lasso_180426"
)

mega_auc <- list()
for (curdir in names(dirSet)) {
	if (curdir %in% c("lasso","pamr")) rngMax <- 20
	else if (curdir %in% "prune") rngMax <- 14
	else rngMax <- 15

	cat(sprintf("***** %s *****\n", curdir))
	dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
	settypes <- c("clinical","mir","rna","cnv","dnam",
		"clinicalArna","clinicalAmir","clinicalAdnam","clinicalAcnv","all")
	ctr <- 1
	auc_set <- list()
	for (settype in settypes) {
	rngDir <- sprintf("%s/rng%i", dataDir,1:rngMax)
	cat(sprintf("Got %i rng files\n",length(rngDir)))
	
		cutoff <- 9
		c7 <- sprintf("%s/%s/predictionResults.txt",
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
		postscript("tmp.eps")
		x <- plotPerf(c7,c("SURVIVEYES","SURVIVENO"))
		dev.off()
	
		y1 <- unlist(lapply(x,function(i) i$auroc))
		auc_set[[settype]] <- y1
	}
	mega_auc[[curdir]] <- unlist(lapply(auc_set,mean))
}

dt <- format(Sys.Date(),"%y%m%d")
pdf(sprintf("OV_%s.pdf",dt)); boxplot(mega_auc,las=1,cex.axis=1.8,main="OV",cex.main=2); dev.off()
