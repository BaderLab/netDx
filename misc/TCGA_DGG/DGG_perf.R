# assess performance of BRCA classifier for different combinations of
# input data types
rm(list=ls())
require(netDx)

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_DGG/output/integrate_170118"

combSet <- c("xpr","dnam","all")	
cols <- c(rep("red",2),"orange")

cat(sprintf("Got %i combs\n", length(combSet)))

out <- list()
overall_acc <- numeric()
for (cur in combSet) {
	inf <- sprintf("%s/%s/predictionResults.txt",inDir,cur)
	dat <- read.delim(inf,sep="\t",h=T,as.is=T)
	out[[cur]] <- perfCalc_multiClass(dat$STATUS,dat$PRED_CLASS)*100
	overall_acc <- c(overall_acc, 
					 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)
}
names(overall_acc) <- combSet

f1 <- unlist(lapply(out, function(x) x[nrow(x),7]))
acc <- unlist(lapply(out, function(x) x[nrow(x),8]))
ppv <- unlist(lapply(out, function(x) x[nrow(x),5]))
par(mfrow=c(2,3),mar=c(9,6,2,2),las=2,cex.axis=1.3,cex.lab=1.7,cex.main=2)
barplot(f1,main="mean F1",col=cols,ylab="F1",ylim=c(0,100))
barplot(acc,main="mean accuracy",col=cols,ylab="accuracy",ylim=c(0,100))
barplot(ppv,main="mean PPV",col=cols,ylab="PPV",ylim=c(0,100))
barplot(overall_acc,main="Accuracy",col=cols,ylab="accuracy",ylim=c(0,100))

