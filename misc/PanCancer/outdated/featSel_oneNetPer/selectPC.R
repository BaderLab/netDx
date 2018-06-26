

	xprFile <- "/home/shraddhapai/BaderLab/2017_PanCancer/GBM/input/GBM_mRNA_core.txt"
xpr <- read.delim(xprFile,sep="\t",h=T,as.is=T)
sname <- xpr[,1]; xpr<- xpr[,-1]

xpr <- t(xpr)
xpr <- xpr[-nrow(xpr),]
class(xpr) <- "numeric"
require(pcaMethods)

