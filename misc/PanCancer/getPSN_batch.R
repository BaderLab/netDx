
curd <- getwd()
rm(list=ls())
cat("multiple clinical nets\n")
setwd("featSel_pathways/pathways")
source("pathways_getPSN.R")
setwd(curd)
rm(list=ls())

cat("multiple clinical nets\n")
setwd("../clinNets")
source("getPSN.R")
setwd(curd)
rm(list=ls())

cat("one clinical net\n")
setwd("../../featSel_oneNetPer")
source("KIRC_clin_oneNetPer_PSN.R")
setwd(curd)

datRoot <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival"
inFile <- list(
	pathways="pathwaysOnly_170502",
	clinNets="clinNets_170430",
	oneClinNet="oneNetPer_FeatSel"
)

# compile Dijkstra
tb <- NULL
for (k in 1:length(inFile)) {
	dat <- read.delim(sprintf("%s/%s/%s_MEAN_Dijkstra_PSN.txt",
		datRoot,inFile[[k]],names(inFile)[k]),sep="\t",h=T,as.is=T)
	if (k==1) {
		tb <- dat 
		colnames(tb) <- colnames(dat)
	} else {
		tb <- rbind(tb, dat)
	}
}
rownames(tb) <- names(inFile)
write.table(tb,file=sprintf("%s/PSN_batch.txt",datRoot),sep="\t",
	col=T,row=T,quote=F)
