#' plot GBM results with multiple CV cutoffs
require(netDx)
require(reshape2)

LUSC_getRes <- function() {
mainD <-  "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output"
dirSet <- list(
#	base="noPrune_180423",
#	baserep="noPrune_sp0.3_180511",
#	baserep_new="noPrune_sp0.3_180527",
	baserep_new2="noPrune_sp0.3_maxEdge6_maxInt40_180529",
#	plasso_303k="pearscale_lasso_topClin1_max30_3K_180528",
#	plasso_306k="pearscale_lasso_topClin1_max30_180528",
#	lasso="lasso_180426",
#	lassoGenes="lassoGenes_180426",
#	pamrGenes_sp2="pamrGenes_180427",
#	pamrGenes_sp1="pamrGenes_sp1_180427",
#	rbfclean="rbfclean_0.20_180507",
	#euc6K="eucclean_180504",
#	eucimpute="eucscale_impute_180507",
#	pearscale="pearscale_180507",
	
#plassoc1="pearscale_lasso_topClin1_180509",
#	plassoc1_new="pearscale_lasso_topClin1_180528",
	lassoc1_new="pearscale_lasso_topClin1_max40_6K_180529"
)
settypes <- c("clinical","mir","rna","prot","cnv",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAcnv","all")

mega_auc <- list()
auc_var <- list()
for (curdir in names(dirSet)) {
dataDir <- sprintf("%s/%s",mainD,dirSet[[curdir]])
#	if (curdir %in% "base") rngMax <- 15
#	else if (curdir %in% "baserep") rngMax <- 20
 rngMax <- 20

auc_set <- list()
for (settype in settypes) {
	if (is.null(auc_var[[settype]])) { auc_var[[settype]] <- c() }

	if (curdir %in% "lassoGenes") {
	rngDir <- paste(sprintf("%s/rng",dataDir), 3:rngMax,sep="")
	#} else if (any(grep("euc",curdir))) {
	#rngDir <- paste(sprintf("%s/rng",dataDir), 1:9,sep="")
	} else {
	rngDir <- paste(sprintf("%s/rng",dataDir), 1:rngMax,sep="")
	}


for (cutoff in 9) {
	if (curdir %in% c("baserep","baserep_new","baserep_new2")) {
		c7 <- sprintf("%s/%s/predictionResults.txt",
				  rngDir,settype)
	} else {
		c7 <- sprintf("%s/%s/cutoff%i/predictionResults.txt",
				  rngDir,settype,cutoff)
	}
	torm <- c()
	for (idx in 1:length(c7)) {
		if (file.exists(c7[idx])) {
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
auc_var[[settype]] <- c(auc_var[[settype]], sd(unlist(auc_set))/sqrt(length(auc_set)))
}
dt <- format(Sys.Date(),"%y%m%d")
pdf(sprintf("LUSC_%s.pdf",dt),width=24,height=6); 
boxplot(mega_auc,main="LUSC",cex.axis=1.7,cex.main=2,las=1); 
abline(h=0.5)
dev.off()

return(mega_auc)
}
