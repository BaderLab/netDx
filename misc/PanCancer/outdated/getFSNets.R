
inDir <- "/Users/shraddhapai/Dropbox (Bader Lab)/netdx/baderlab/2017_pancancer_survival/pathwaysOnly_170502/featSelNets"
inFiles <- c("pathways_thresh10_pctPass0.70_SURVIVENO_netScores.txt",
		"pathways_thresh10_pctPass0.70_SURVIVEYES_netScores.txt")

for (f in inFiles) {
  dat <- read.delim(sprintf("%s/%s",inDir,f),sep="\t",h=T,as.is=T)
	nets<- dat[,1]; dat <- dat[,-1]
	dat <- dat >= 10
	idx <- rowSums(dat,na.rm=TRUE)>0.7*ncol(dat)
	
	cat(sprintf("%s\n", f))
	print(nets[idx])	
	cat("------\n")
}
