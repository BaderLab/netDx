#' plot LUSC results with multiple CV cutoffs
rm(list=ls())

#dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output/multiCutoff_180119"
#dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/prune_180125"
#dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/pruneCheckIntegr_180126"
dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/pruneRBF_180130"
#dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/prunePCA_180126"

rngDir <- paste(sprintf("%s/rng",dataDir), 1:10,sep="")

require(netDx)
for (cutoff in 7:9) {
	c7 <- sprintf("%s/clinical/cutoff%i/predictionResults.txt",rngDir,cutoff)
	torm <- c()
	for (idx in 1:length(c7)) {
		dat <- read.delim(c7[idx],sep="\t",h=T,as.is=T)
		x1 <- sum(dat$STATUS=="SURVIVEYES")
		x2 <- sum(dat$STATUS=="SURVIVENO")
	#	cat(sprintf("%i: %i YES, %i NO\n", idx,x1,x2))
		if (x1<1 & x2<1) torm <- c(torm, idx)
	}
	cat(sprintf("%i: removing %i\n", cutoff,length(torm)))
	if (length(torm)>0) c7 <- c7[-torm]
	postscript(sprintf("LUSC_cutoff%i.eps",cutoff));
	x <- plotPerf(c7,c("SURVIVEYES","SURVIVENO"))
	y <- unlist(lapply(x,function(i) i$auroc))
	cat(sprintf("%i, mean auroc= %1.2f\n", cutoff, mean(y)))
	dev.off()
}

###
