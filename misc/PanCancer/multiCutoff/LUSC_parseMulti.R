#' plot GBM results with multiple CV cutoffs
rm(list=ls())
require(netDx)
require(reshape2)

dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/LUSC/output/prunedPearson_180212"

settypes <- c("clinical","mir","rna","prot","cnv",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAcnv","all")
outmat <- matrix(NA,nrow=length(settypes),ncol=9)
meas <- paste(rep(7:9,each=3),c("auroc","aupr","accuracy"),sep="_")
rownames(outmat)<- settypes
colnames(outmat) <- meas
ctr <- 1
outD <- sprintf("LUSC_%s",basename(dataDir))
if (!file.exists(outD)) dir.create(outD)

auc_set <- list()
for (settype in settypes) {
###	if (settype %in% "clinicalArna") 
###		dataDir <- dataDir_both
###	else 
###		dataDir <- dataDir_each
	rngDir <- paste(sprintf("%s/rng",dataDir), 1:93,sep="")

colctr <- 1
for (cutoff in 9) {
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
	postscript(sprintf("%s/%s_cutoff%i.eps",outD,settype,cutoff));
	x <- plotPerf(c7,c("SURVIVEYES","SURVIVENO"))
	dev.off()

	y1 <- unlist(lapply(x,function(i) i$auroc))
	y2 <- unlist(lapply(x,function(i) i$aupr))
	y3 <- unlist(lapply(x,function(i) i$accuracy))
	outmat[ctr,colctr+(0:2)] <- c(mean(y1),mean(y2),mean(y3))
	auc_set[[settype]] <- y1

	colctr <- colctr+3
}
ctr <- ctr+1
}
print(round(outmat,digits=2))

pdf("lusc_auc.pdf",width=13,height=5);
 boxplot(auc_set,cex.axis=0.6,pars=list(boxwex=0.3)); 
	abline(h=median(auc_set[["clinical"]]));
 barplot(unlist(lapply(auc_set,mean)),las=1,cex.axis=1.3,font.axis=2,
	ylim=c(0.5,1),main="LUSC")
dev.off()

write.table(round(outmat,digits=2),file=sprintf("%s/perf.txt",outD),sep="\t",
			col=T,row=T,quote=F)


