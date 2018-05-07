#' plot GBM results with multiple CV cutoffs
require(netDx)
require(reshape2)

#dataDir <- "/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pruneTrain_180419"

KIRC_getRes <- function() {
mainD <-  "/home/shraddhapai/BaderLab/PanCancer_KIRC/output"

dirSet <- list(
	base="noPrune_180423",
	lasso="lasso_180426",
	pamr="pamr_180426",
	euc6K="eucclean_180503"
#	ridge="ridgeAbsFix_180426"
)

settypes <- c("clinical","mir","rna","prot","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAdnam",
	"clinicalAcnv","all")
mega_auc <- list()

for (curdir in names(dirSet)) {
dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
	rngMax <- 20
	if (any(grep("base",curdir))) {
		rngMax <- 15
	}	

auc_set <- list()
for (settype in settypes) {
	rngDir <- paste(sprintf("%s/rng",dataDir), 1:rngMax,sep="")

colctr <- 1
	cutoff <- 9
	if (curdir=="euc6K") {
		c7 <- sprintf("%s/%s/cutoff9/predictionResults.txt",
				  rngDir,settype,cutoff)
	} else {
		c7 <- sprintf("%s/%s/predictionResults.txt",
				  rngDir,settype,cutoff)
	}
	
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
	}
	mega_auc[[curdir]] <- unlist(lapply(auc_set,mean))
}
pdf(sprintf("KIRC_%s.pdf",format(Sys.Date(),"%y%m%d"))); 
boxplot(mega_auc,main="KIRC",cex.axis=1.7,cex.main=2,las=1); dev.off()

return(mega_auc)
}


