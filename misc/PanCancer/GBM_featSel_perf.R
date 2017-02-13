# assess performance of BRCA classifier for different combinations of
# input data types
rm(list=ls())
require(netDx)

#inDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_noMut_170208"
inDir <- "/mnt/data2/BaderLab/PanCancer_GBM/output/featSel_170209"

cols <- c(brewer.pal(n=4,name="Dark2"),"red")

NBT_val <- 0.705 # best GBM AUCROC from Yuan et al NBT paper Figure 4.

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

mega <- list()
for (rngSeed in 1:30) {
	out <- list()
	overall_acc <- numeric()
	curRoc	<- list()
	curPr	<- list()
	cur <- "base"

	# collect info for ROC curve
	inf <- sprintf("%s/rng%i/predictionResults.txt",
				   inDir,rngSeed,cur)
	dat <- read.delim(inf,sep="\t",h=T,as.is=T)
	pred <- prediction(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE,
					  dat$STATUS=="SURVIVEYES")

	c1 <- "SURVIVEYES" #numc[1]
	tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
	tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
	fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
	fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)

	# dummy score column needed for perfCalc()
	tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn) 
	curRoc[[cur]] <- performance(pred,"tpr","fpr")
	curPr[[cur]] <- performance(pred,"prec","rec")

	out[[cur]] <- perfCalc(tmp)
	# add AUC-PR and AUCROC
	out[[cur]]$auc <- performance(pred, "auc")@y.values[[1]]
	out[[cur]]$prauc <- prauc(curPr[[cur]])
	overall_acc <- c(overall_acc, 
		 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)
	names(overall_acc) <- "base"
	mega[[as.character(rngSeed)]] <- list(out=out,overall_acc=overall_acc,
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

# extract performance measures, compile to plotting format
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
pdf("GBM_perf.pdf",width=8,height=4)
tryCatch({

par(las=2,cex.axis=1.1,cex.lab=1.5,cex.main=1.3)
stats <- list(F1=f1,PPV=ppv,Accuracy=acc,AUCROC=aucroc,AUCPR=aucpr)

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
barx <- x
#text(x=colMeans(x),y=0.9, sprintf("p < %1.0e",pvalues))

for (i in 1) {
	sds <- rbind(sd(ppv[,i]),sd(acc[,i]),
		 sd(aucroc[,i]), sd(aucpr[,i]))
	sds <- sds/sqrt(nrow(ppv))
	segments(x0=x[i,], x1=x[i,],y0=mu[,i]-sds,y1=mu[,i]+sds)
}
#plot(0,0,bty='n',xaxt='n',yaxt='n',type='n',xlab="",ylab="")
legend("topright", legend=colnames(ppv),fill=cols,bty='n')

cat("Best performance:\n")
out <- list()
for (nm in names(stats)) {
	cur <- stats[[nm]];
	x <- unlist(sapply(1:ncol(cur), function(x) max(cur[,x])))
	names(x) <- colnames(cur)
	out[[nm]] <- x
}
maxstat <- do.call("rbind",out)
print(signif(maxstat,2))
# plot max values with 'x's
points(x=barx,y=maxstat[-which(rownames(maxstat) %in% "F1")],
	pch=4,cex=1.5,col='red')
# mark NBT values 
points(x=barx[3],y=NBT_val,pch=4,col='blue',cex=1.5)

}, error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

cat("Average performance:\n")
out <- list()
for (nm in names(stats)) {
	cur <- stats[[nm]];
	x <- unlist(sapply(1:ncol(cur), function(x) mean(cur[,x])))
	names(x) <- colnames(cur)
	out[[nm]] <- x
}
mustat <- do.call("rbind",out)
print(signif(mustat,2))


# plot individual ROC curves
pdf("GBM_perf_ROC.pdf",width=11,height=4)
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


pullROC <- function(str) {
	z <- lapply(mega, function(x) { y <- x$roc; y[[str]]})
	z
}

pdf("GBM_aggROC.pdf",width=10,height=5)
source("plotROC_multi.R"); 
par(mfrow=c(2,3))
for (str in colnames(maxstat)) {
	x <- pullROC(str)
	plotROC_multi(x,which.max(stats[["AUCROC"]][,str]))
	title(str)
}
dev.off()

