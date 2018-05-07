#' plot GBM results for kernel variations

require(netDx)
require(reshape2)

GBM_getRes <- function() {
mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output"
dirSet <- list(
	base="noPrune_180423",
	ridge_fix="ridge_AbsFix_180426",
	lassoGenes_sp1="lassoGenes_incClin_180426",
	pamrGenes="pamrGenes_incClin_180427",
	#rbf0.05="lassoUni_rbf_0.05",
	rbf0.1="lassoUni_rbf_0.1_180502",
	rbf0.25="lassoUni_rbf_0.25_180502",
	#euc_1K="eucscale_sp2max1000_180503",
	euc_6K="eucclean_180503",
	euc_6K_group="eucscale_sp2max6000_grouped_180503"
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

	if (any(c(grep("lasso",curdir),grep("ridge",curdir)))) {
		rngDir <- paste("rng",1:18,sep="")
	} else if (any(c(grep("rbf0.1",curdir)))){
		rngDir <- paste("rng",1:8,sep="")
	} else if (any(c(grep("rbf0.25",curdir)))){
		rngDir <- paste("rng",1:8,sep="")
	} else if (any(c(grep("euc_1K",curdir)))){
		rngDir <- paste("rng",1:12,sep="")
	} else if (curdir =="euc_6K"){
		rngDir <- paste("rng",1:20,sep="")
	} else if (curdir =="euc_6K_group"){
		rngDir <- paste("rng",1:20,sep="")
	} else {
	rngDir <- dir(path=dataDir,pattern="rng")
	}
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
pdf(sprintf("GBM_%s.pdf",dt),width=18,height=6);
boxplot( mega_auc,las=1,cex.axis=1.7,cex.main=2,main="GBM",
	at=1:length(mega_auc)); 
tmp <- unlist(numSplits)
text(1:length(mega_auc),0.5,sprintf("N=%i",tmp))
abline(h=0.5)
dev.off()

return(mega_auc)
}
