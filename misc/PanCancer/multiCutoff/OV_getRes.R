#' plot GBM results with multiple CV cutoffs
require(netDx)
require(reshape2)

OV_getRes <- function() {
mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/OV/output"

maxRng <- 100
settypes <- c("clinical","mir","rna","prot","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAdnam",
	"clinicalAcnv","all")
dirSet <- list(
#	base="noPrune_180423",
	#baserep="noprune_sp0.3_180511",
	baserep="noprune_sp0.3_180527",
	baserep1="noprune_sp1_180512",
#	prune="pruneTrain_180419",
#	lasso="lasso_180426",
#	euc6K="eucscale_180504",
	euc_correct="eucscale_180528"
#	rbfclean="rbfclean_0.2_180507"
)

mega_auc <- list()
for (curdir in names(dirSet)) {
#	if (curdir %in% c("lasso","pamr","euc6K","rbfclean")) rngMax <- 20
#	else if (curdir %in% "prune") rngMax <- 14
	
#	if (curdir=="base") rngMax <- 15 
#	else if (curdir=="baserep") rngMax <- 20
	rngMax <- 20

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
		if (curdir %in% c("euc6K","rbfclean","euc_correct")) {
		c7 <- sprintf("%s/%s/cutoff9/predictionResults.txt",
					  rngDir,settype,cutoff)
		} else {
		c7 <- sprintf("%s/%s/predictionResults.txt",
					  rngDir,settype,cutoff)
		}
		torm <- c()
		for (idx in 1:length(c7)) {
			if (file.exists(c7[idx])){
			dat <- read.delim(c7[idx],sep="\t",h=T,as.is=T)
			x1 <- sum(dat$STATUS=="SURVIVEYES")
			x2 <- sum(dat$STATUS=="SURVIVENO")
			if (x1<1 & x2<1) torm <- c(torm, idx)
			} else {
				torm <- c(torm,idx)
			}
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
pdf(sprintf("OV_%s.pdf",dt),width=12,height=6); 
boxplot(mega_auc,las=1,cex.axis=1.8,main="OV",cex.main=2); dev.off()

	return(mega_auc)

}
