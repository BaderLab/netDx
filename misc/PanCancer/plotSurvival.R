#' plot survival boxplot

inRoot <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/"

inSet <- c("KIRC","LUSC","OV","GBM")

pCutoff <- list(KIRC=4,LUSC=2,OV=3,GBM=1)

require(ggplot2)
dat <- list()
pList <- list()
for (k in inSet) {
	surv <- read.delim(sprintf("%s/2017_TCGA_%s/input/%s_binary_survival.txt",
		inRoot,k,k),sep="\t",h=T,as.is=T)
	pheno <- read.delim(sprintf("%s/2017_TCGA_%s/input/%s_OS_core.txt",
		inRoot,k,k),sep="\t",h=T,as.is=T)
	x <- merge(x=surv,y=pheno,by="feature")
	x$OS_years <- x$OS_OS/365
	x$OS_mos <- x$OS_OS/30
	p <- ggplot(x,aes(factor(is_alive),OS_years)) + geom_boxplot()
	p <- p + geom_hline(yintercept=pCutoff[[k]],col='red')
	p2 <- wilcox.test(x$OS_years[which(x$is_alive<1)],
										x$OS_years[which(x$is_alive>0)],"less")$p.value
	p <- p + ggtitle(sprintf("%s: p < %1.2e",k,p2))

	mu1 <- median(x$OS_mos[which(x$is_alive<1)])
	mu2 <- median(x$OS_mos[which(x$is_alive>0)])
	cat(sprintf("%s: Low = %1.1f mos; High = %1.1f mos; diff= %1.1f YEARS\n",
		k,mu1,mu2,(mu2-mu1)/12))
	pList[[k]] <- p
}

source("multiplot.R")
multiplot(plotlist=pList,layout=matrix(1:4,ncol=2))
