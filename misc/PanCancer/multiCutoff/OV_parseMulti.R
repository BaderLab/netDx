#' plot GBM results with multiple CV cutoffs
rm(list=ls())
require(netDx)
require(reshape2)

#dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/OV/output/prune_180423"
dataDir <- "/home/shraddhapai/BaderLab/2017_PanCancer/OV/output/ridge_180420"

maxRng <- 100
settypes <- c("clinical","mir","rna","prot","cnv","dnam",
	"clinicalArna","clinicalAmir","clinicalAprot","clinicalAdnam",
	"clinicalAcnv","all")
outmat <- matrix(NA,nrow=length(settypes),ncol=3)
meas <- paste(rep(9,each=3),c("auroc","aupr","accuracy"),sep="_")
rownames(outmat)<- settypes
colnames(outmat) <- meas
ctr <- 1
outD <- sprintf("OV_%s",basename(dataDir))
if (!file.exists(outD)) dir.create(outD)

var_set <- list()
auc_set <- list()
for (settype in settypes) {
###	if (settype %in% "clinicalArna") 
###		dataDir <- dataDir_both
###	else 
###		dataDir <- dataDir_each
rngDir <- paste(sprintf("%s/rng",dataDir), 1:maxRng,sep="")

colctr <- 1
for (cutoff in 9) {
	c7 <- sprintf("%s/%s/predictionResults.txt",
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
	tmp <- c()
	cur <- auc_set[[settype]]
	for (k in 3:length(cur)) {
		tmp <- c(tmp, sd(cur[1:k]))
	}
	var_set[[settype]] <- data.frame(type=settype,numsplits=4:length(cur),
			pctChangeVar=diff(tmp^2)/(tmp[-1]^2))
}
ctr <- ctr+1
}

#auc_set <- auc_set[which(names(auc_set)%in% c("clinical","rna","clinicalArna","all"))]
cat("-----------------\n")
cat(sprintf("OV: Base dir: %s\n", basename(dataDir)))
cat(sprintf("%i splits\n", length(auc_set[[1]])))
cat("-----------------\n")
print(round(outmat,digits=2))
cat("-----------------\n")
pdf("ov_auc.pdf",width=16,height=5);
 boxplot(auc_set,cex.axis=0.6,pars=list(boxwex=0.3)); 
	abline(h=median(auc_set[["clinical"]]));

mu <- unlist(lapply(auc_set,mean))
err <- unlist(lapply(auc_set,sd))
xpos <-  barplot(mu,las=1,cex=1.3,cex.names=0.8,
	ylim=c(0,1),main="OV")
segments(x0=xpos,y0=mu-err,y1=mu+err)

wmw <- wilcox.test(auc_set[["all"]],auc_set[["clinicalArna"]],alternative="greater")
cat(sprintf("All > clin+rna: p < %1.2e\n", wmw$p.value))

dev.off()

write.table(round(outmat,digits=2),file=sprintf("%s/perf.txt",outD),sep="\t",
			col=T,row=T,quote=F)

# plot SEM as function of num rounds
setName <- sprintf("OV_%s", basename(dataDir))
for (settype in settypes) {
	tmp <- var_set[[settype]][,3]
	tmp2 <- c()
	for (m in 1:(length(tmp)-2)) tmp2 <- c(tmp2,mean(tmp[m:(m+2)]))
	cat(sprintf("%s: < 1%% change: %i\n", settype,min(which(abs(tmp2) < 0.01))))
}
var_set <- do.call("rbind",var_set)
require(ggplot2)
p <- ggplot(var_set,aes(x=numsplits,y=pctChangeVar)) 
p <- p+ geom_smooth(aes(colour=type),method="loess",span=0.1,se=FALSE,lwd=0.5,
	alpha=0.5)
p <- p + ggtitle(setName) + ylim(c(-0.25,0.25))
p <- p + geom_vline(xintercept=c(10,15,25),lty=3)
pdf(sprintf("%s.pdf",setName),width=8,height=3); print(p);dev.off()

