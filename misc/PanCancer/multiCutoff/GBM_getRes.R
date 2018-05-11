#' plot GBM results for kernel variations

require(netDx)
require(reshape2)

GBM_getRes <- function() {
mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output"
dirSet <- list(
	base="noPrune_180423",
#	ridge_fix="ridge_AbsFix_180426",
#	lassoGenes_sp1="lassoGenes_incClin_180426",
#	pamrGenes="pamrGenes_incClin_180427",
#	euc_6K="eucclean_180503",
	eucimpute="eucclean_impute_180507",
	#rbfclean="rbfclean_0.20_180507",
#	pearscale="pearscale_180507",
	pearimpute="pearscale_impute_180508"
#	pimp_40c2="pearimp_topX40_topClin2_180509",
#	pimp_40c3="pearimp_topX40_topClin3_180509",
#	pimp_100c3="pearimp_topX100_topClin3_180509",
#	pimp_200c3="pearimp_topX200_topClin3_180509",
#	pimp_20c3="pearimp_topX20_topClin3_180509",
#	pimp_30c3="pearimp_topX30_topClin3_180509"
#	eimp_100c3="eucimp_topX100_topClin3_180509",
#	eimp_50c3="eucimp_topX50_topClin3_180509",
#	eimp_20c3="eucimp_topX20_topClin3_180509"
	#rbf5="lassoUni_rbf_5",
	#rbf10="lassoUni_rbf_10"
)

mega_auc <- list()
numSplits <- list()
for (curdir in names(dirSet)) {
cat(sprintf("***** %s *****\n", curdir))
dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
settypes <- c("clinical","mir","rna","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAdnam","clinicalAcnv","all")
ctr <- 1

auc_set <- list()
for (settype in settypes) {
	#print(dataDir)
cutoff <-9

#	if (any(c(grep("lasso",curdir),grep("ridge",curdir)))) {
#		rngDir <- paste("rng",1:18,sep="")
#	} else if (any(c(grep("rbf0.1",curdir)))){
#		rngDir <- paste("rng",1:8,sep="")
#	} else if (any(c(grep("rbf0.25",curdir)))){
#		rngDir <- paste("rng",1:8,sep="")
#	} else if (any(c(grep("euc_1K",curdir)))){
#		rngDir <- paste("rng",1:12,sep="")
#	} else if (curdir =="euc_6K"){
#		rngDir <- paste("rng",1:20,sep="")
#	} else if (any(grep("pimp",curdir))){
#		rngDir <- paste("rng",1:14,sep="")
#	} else {
	if (curdir=="base") rngMax<- 15 else rngMax <- 20
	rngDir <- paste("rng",1:rngMax,sep="") #dir(path=dataDir,pattern="rng")
#	}
	numSplits[[curdir]] <- length(rngDir)

	cat(sprintf("Got %i rng files\n",length(rngDir)))
	rngDir <- sprintf("%s/%s",dataDir,rngDir)
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

	postscript("tmp.eps")
	x <- plotPerf(c7,c("SURVIVEYES","SURVIVENO"))
	dev.off()

	y1 <- unlist(lapply(x,function(i) i$auroc))
	auc_set[[settype]] <- y1
ctr <- ctr+1
}
mega_auc[[curdir]] <- unlist(lapply(auc_set,mean))

}

dt <- format(Sys.Date(),"%y%m%d")
require(gplots)
pdf(sprintf("GBM_%s.pdf",dt),width=24,height=6);
boxplot( mega_auc,las=1,cex.axis=1.7,cex.main=2,main="GBM",
	at=1:length(mega_auc),cex.lab=0.8); 
tmp <- unlist(numSplits)
text(1:length(mega_auc),0.5,sprintf("N=%i",tmp))
abline(h=0.5)
dev.off()

return(mega_auc)
}
