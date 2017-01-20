# assess performance of BRCA classifier for different combinations of
# input data types
rm(list=ls())
require(netDx)

#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/output/run170119"
inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/output/run170119_noSMC_MCI"

combSet <- c("clinical","genetic","imaging","all")
cols <- c(rep("red",3),"orange")

cat(sprintf("Got %i combs\n", length(combSet)))

mega <- list()
for (rngSeed in 1:10) {
	out <- list()
	overall_acc <- numeric()
	for (cur in combSet) {
		inf <- sprintf("%s/rng%i/%s/predictionResults.txt",
					   inDir,rngSeed,cur)
		dat <- read.delim(inf,sep="\t",h=T,as.is=T)
		c1 <- "worse" #numc[1]
		tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
		tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
		fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
		fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)

		# dummy score column needed for perfCalc()
		tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn) 
		out[[cur]] <- perfCalc(tmp)
		overall_acc <- c(overall_acc, 
						 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)
	}
	names(overall_acc) <- combSet
	mega[[rngSeed]] <- list(out=out,overall_acc=overall_acc)
}

# get stat interest from data structure with all results
func <- function(str) {
	tmp <- lapply(mega, function(x) {
		y <- x$out
		z <- unlist(lapply(y,function(w) { w$stats[[str]]}))
		z
	})
	tmp <- do.call("rbind",tmp)

	# normalize to first column
	#for (k in ncol(tmp):1) tmp[,k] <- tmp[,k]-tmp[,1]
	tmp
}

f1 <- func("f1")
ppv <- func("ppv")
acc <- lapply(mega, function(x) { x$overall_acc/100 })
acc <- do.call("rbind",acc)
#for (k in ncol(acc):1) acc[,k] <- acc[,k]-acc[,1]

# plot
par(mfrow=c(2,2),mar=c(6,4,2,2),las=2,
	cex.axis=1.1,cex.lab=1.5,cex.main=1.3)
barplot(f1,beside=TRUE, main="F1",ylab="")
barplot(ppv,beside=TRUE, main="PPV",ylab="")
barplot(acc,beside=TRUE,main="Accuracy",ylab="")

mu <- rbind(colMeans(f1), colMeans(ppv), colMeans(acc))
sds <- rbind(sd(f1[,4]), sd(ppv[,4]),sd(acc[,4]))
rownames(mu) <- c("F1","PPV","ACC")
x <- barplot(t(mu),col=c(rep("red",3),"orange"),beside=TRUE,
		#main="Mean over 10 runs\n(rel. to clinical only)",ylab="delta clinical",
		main="Mean over 10 runs",ylab="",
		#ylim=c(-.1,+0.05),cex.main=1.1,cex.lab=1.2)
		cex.main=1.1,cex.lab=1.2,ylim=c(0.4,0.7))
segments(x0=x[4,], x1=x[4,],y0=mu[,4]-sds,y1=mu[,4]+sds)

