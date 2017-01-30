# assess performance of BRCA classifier for different levels of missing
# data
rm(list=ls())
require(netDx)
require(RColorBrewer)


inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/TCGA_breastCancer/output/missing_170118"

combSet <- paste("miss",c(10,50,70,85,90,95,99),sep="")
cols <- brewer.pal(n=length(combSet),name="Blues")

cat(sprintf("Got %i combs\n", length(combSet)))

out <- list()
overall_acc <- numeric()
pctMiss <- numeric()
numInt <- numeric()
for (cur in combSet) {
	inf <- sprintf("%s/%s/predictionResults.txt",inDir,cur)
	dat <- read.delim(inf,sep="\t",h=T,as.is=T)
	dat <- dat[-which(dat$STATUS %in% "Normal"),]
	out[[cur]] <- perfCalc_multiClass(dat$STATUS,dat$PRED_CLASS)*100
	overall_acc <- c(overall_acc, 
					 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)

	tmp <- system(sprintf("wc -l %s/%s/tmp/INTERACTIONS/1.1.txt",
							 inDir,cur),intern=TRUE)
	tmp <- strsplit(tmp," ")[[1]][4]
	numInt <- c(numInt,as.integer(tmp))

	dat <- read.delim(sprintf("%s/%s/profiles/xpr.profile",inDir,cur),
					  h=T,as.is=TRUE)
	pctMiss <- c(pctMiss, sum(is.na(dat))/(nrow(dat)*ncol(dat)))
}
names(overall_acc) <- combSet

tot <- unlist(lapply(out,function(x) sum(x[1,1:4])/100))
f1 <- unlist(lapply(out, function(x) x[nrow(x),7]))
acc <- unlist(lapply(out, function(x) x[nrow(x),8]))
ppv <- unlist(lapply(out, function(x) x[nrow(x),5]))

par(mfrow=c(3,3),mar=c(9,6,2,2),las=2,cex.axis=1.3,cex.lab=1.7,cex.main=2)
barplot(f1,main="mean F1",col=cols,ylab="F1",ylim=c(0,100))
barplot(acc,main="mean accuracy",col=cols,ylab="accuracy",ylim=c(0,100))
abline(h=c(50,85),col='red',lwd=2,lty=2)
barplot(ppv,main="mean PPV",col=cols,ylab="PPV",ylim=c(0,100))
abline(h=50,col='red',lwd=2,lty=2)
barplot(overall_acc,main="Accuracy",col=cols,ylab="accuracy",ylim=c(0,100))
abline(h=c(50,85),col='red',lwd=2,lty=2)
barplot(tot,main="total classified")
names(numInt) <- names(tot)
barplot(numInt/1000,main="# interactions (x1000)")
names(pctMiss) <- names(tot)
barplot(pctMiss*100,main="% missing in profile",ylim=c(0,100))

