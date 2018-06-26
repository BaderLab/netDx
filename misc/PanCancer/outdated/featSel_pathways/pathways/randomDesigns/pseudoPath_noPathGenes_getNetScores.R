#' count num times nets come up in randomly sampled set

inDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/randomD_pseudoPath_noPathGenes_170804"

rngDir <- dir(inDir,pattern="rng")

for (gp in c("SURVIVEYES","SURVIVENO")) {
	netSet <- list()
	for (cur in rngDir) {
		netSet[[cur]] <- read.delim(sprintf("%s/%s/%s/netNames.txt",
			inDir,cur,gp),sep="\t",h=F,as.is=T)[,1]
	}
	netSet <- unlist(netSet)
	netSc <- table(netSet)
	netSc <- data.frame(netName=names(netSc), netScore=as.integer(netSc))
	netSc <- netSc[order(netSc[,2],decreasing=TRUE),]
	print(head(netSc))
	write.table(netSc,file=sprintf("%s/%s.netTally.txt",inDir,gp),
		sep="\t",col=T,row=F,quote=F)
}
