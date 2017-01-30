# assess performance of BRCA classifier for different combinations of
# input data types
rm(list=ls())
require(netDx)

#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/output/run170119"
inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_ADNI/output/featSel_170125_test2"

combSet <- "all" #c("clinical","genetic","imaging","all")
cols <- c(brewer.pal(n=3,name="Dark2"),"red")

# --------------------------------------------------------------
# helper routines

require(RColorBrewer)
.plotROC <- function(dat,pType) {
	plot(NA,NA,
		 ylab="",xlab="",xlim=c(0,1),ylim=c(0,1),bty='n')
	for (k in 1:length(dat)) {
		if (k==4) lwd=1.4 else lwd=1
		x <- dat[[k]]@x.values[[1]]
		y <- dat[[k]]@y.values[[1]]
		points(x,y,type='l',col=cols[k],lwd=lwd)
	}
	text(0.5,0,sprintf("N=%i",length(dat[[1]]@x.values[[1]])))

	box(col='grey80')
	if (pType == "roc") abline(0,1,col='grey50')
	else abline(h=0.5,col='grey50')
}

# given output of performance("precall") compute AUC-PR
prauc <- function(dat) {
	x <- dat@x.values[[1]] # recall
	y <- dat@y.values[[1]] # precision

	# remove NAN
	idx <- which(is.nan(y))
	if (any(idx)) { x <- x[-idx]; y <- y[-idx]}

	#x <- c(0,x,1); y <- c(1,y,0) # make sure points go from 0,0 to 1,1

	pracma::trapz(x,y)
}

# --------------------------------------------------------------
# begin work
cat(sprintf("Got %i combs\n", length(combSet)))

mega <- list()
for (rngSeed in 1){  #:10) {
	out <- list()
	overall_acc <- numeric()
	curRoc	<- list()
	curPr	<- list()
	for (cur in combSet) {
		# collect info for ROC curve
		#x <- outRes$same$roc; 
		#x2 <- outRes$worse$roc;
		#curRoc[[cur]] <- list(same=outRes$same$roc, 
		#					  worse=outRes$worse$roc)
		#curPr[[cur]] <- list(same=outRes$same$precall,
	#						 worse=outRes$worse$precall)

		inf <- sprintf("%s/rng%i/%s/predictionResults.txt",
					   inDir,rngSeed,cur)
		dat <- read.delim(inf,sep="\t",h=T,as.is=T)

		
		out[[cur]] <- perfCalc_multiClass(dat$STATUS,dat$PRED_CLASS)*100
		overall_acc <- c(overall_acc, 
				 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)

		browser()
		#out[[cur]] <- perfCalc(tmp)
		# add AUC-PR and AUCROC
		#out[[cur]]$auc <- performance(pred, "auc")@y.values[[1]]
		#out[[cur]]$prauc <- prauc(curPr[[cur]])
		#out[[cur]]$aucroc #<- list(same=outRes$same$auc, 
						#		  worse=outRes$worse$auc,
						#		  comb=mean(c(outRes$same$auc, 
						#					outRes$worse$auc)))
		#out[[cur]]$aucpr <- list(same=prauc(outRes$same$precall),
		#						 worse=prauc(outRes$worse$precall))
		#out[[cur]]$aucpr$comb <- mean(unlist(out[[cur]]$aucpr))
	}
	names(overall_acc) <- combSet
	mega[[rngSeed]] <- list(out=out,overall_acc=overall_acc,
							roc=curRoc,pr=curPr)
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
aucroc	<- lapply(mega,function(x) {
	z <- lapply(x$out, function(y) y$auc)
	unlist(z)
	})
aucroc <- do.call("rbind",aucroc)
aucpr	<- lapply(mega,function(x) {
	z <- lapply(x$out, function(y) y$prauc)
	unlist(z)
	})
aucpr <- do.call("rbind",aucpr)

# plot
pdf("ADNI_perf.pdf",width=8,height=4)
tryCatch({
#par(mfrow=c(2,4),mar=c(6,4,2,2),las=2,
#	cex.axis=1.1,cex.lab=1.5,cex.main=1.3)
layout(matrix(c(1:4,5,5,6,6),ncol=4,byrow=TRUE))
par(las=2,cex.axis=1.1,cex.lab=1.5,cex.main=1.3)
stats <- list(F1=f1,PPV=ppv,Accuracy=acc,AUCROC=aucroc,AUCPR=aucpr)
require(reshape2)
pvalues <- numeric()
for (nm in names(stats)) {
	barplot(stats[[nm]],beside=TRUE,main=nm,ylab="")
	tmp <- melt(stats[[nm]])
	colnames(tmp)[2] <- "datatype"

	# examine effect of datatype
	fit <- summary(aov(value~datatype,data=tmp))
	p <- fit[[1]][["Pr(>F)"]][1]
	pvalues <- c(pvalues,p)
}
names(pvalues) <- names(stats)
print(pvalues)

mu <- rbind(colMeans(ppv), colMeans(acc),
			colMeans(aucroc),colMeans(aucpr))
rownames(mu) <- c("PPV","ACC","AUCROC","AUCPR")
par(las=1)
x <- barplot(t(mu),col=cols,beside=TRUE,
		#main="Mean over 10 runs\n(rel. to clinical only)",
		#ylab="delta clinical",
		main="Mean over 10 runs",ylab="",
		#ylim=c(-.1,+0.05),cex.main=1.1,cex.lab=1.2)
		cex.main=1.1,cex.lab=1.2,ylim=c(0,1))
text(x=colMeans(x),y=0.9, sprintf("p < %1.0e",pvalues))

for (i in 1:4) {
	sds <- rbind(sd(ppv[,i]),sd(acc[,i]),
			 sd(aucroc[,i]), sd(aucpr[,i]))
	segments(x0=x[i,], x1=x[i,],y0=mu[,i]-sds,y1=mu[,i]+sds)
}
plot(0,0,bty='n',xaxt='n',yaxt='n',type='n',xlab="",ylab="")
legend("right", legend=colnames(ppv),fill=cols,bty='n')

# now normalize by first column
stats2 <- lapply(stats, function(x) {
		x2 <- x
		for (k in ncol(x):1) {
			x2[,k] <- x2[,k]-x2[,1]
		}
		x2
})
mu2 <- lapply(stats2, function(x) colMeans(x))
mu2 <- do.call("rbind",mu2)
layout(matrix(c(1,1,2,3,4,5,6,7),ncol=4,byrow=TRUE))
x <- barplot(t(mu2),col=cols,beside=TRUE,
		#main="Mean over 10 runs\n(rel. to clinical only)",
		#ylab="delta clinical",
		main="Mean (normalized to clinical)",ylab="",
		#ylim=c(-.1,+0.05),cex.main=1.1,cex.lab=1.2)
		cex.main=1.1,cex.lab=1.2,ylim=c(-0.1,0.1))
for (i in 1:4) {
	sds <- lapply(stats2,function(x) sd(x[,i]))
	sds <- do.call("rbind",sds)
	segments(x0=x[i,], x1=x[i,],y0=mu2[,i]-sds,y1=mu2[,i]+sds)
}

}, error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

# plot individual ROC curves
pdf("ADNI_perf_ROC.pdf",width=11,height=4)
par(mfrow=c(2,5),mar=c(3,3,1,1),las=1)
tryCatch({
	lapply(mega,function(x) .plotROC(x$roc,"roc"))
	legend("bottomright", bty='n', legend=names(mega[[1]]$roc),
		   fill=cols)
	lapply(mega,function(x) .plotROC(x$pr,"PR"))
	legend("bottomright", bty='n', legend=names(mega[[1]]$roc),
		   fill=cols)

},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})



