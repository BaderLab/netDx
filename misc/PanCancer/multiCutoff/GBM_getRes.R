#' plot GBM results for kernel variations

rm(list=ls())
require(netDx)
require(reshape2)

mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output"
dirSet <- list(
	base="noPrune_180423",
#	rbf3="rbf0.3_noSex_180424",
#	rbf5="rbf0.5_noSex_180424",
#	rbfother="rbf_noSex_180424",
#	tanh="tanh_noSex_180424"
	lasso_sp2="lassoGenes_180426",
	lasso_sp1="lassoGenes_sparse1_180426"
)

mega_auc <- list()
for (curdir in names(dirSet)) {
cat(sprintf("***** %s *****\n", curdir))
dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
settypes <- c("clinical","mir","rna","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAdnam","clinicalAcnv","all")
ctr <- 1

auc_set <- list()
for (settype in settypes) {
	#print(dataDir)
cutoff <-8

	if (any(grep("lasso",curdir))) {
		rngDir <- paste("rng",1:8,sep="")
	} else {
	rngDir <- dir(path=dataDir,pattern="rng")
	}

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

pdf("test.pdf"); boxplot(mega_auc); dev.off()
