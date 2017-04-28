# plot AUCROC for predictor

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneNetPer_FeatSel"

inDir <- sprintf("%s/fromAhmad_170428/all_rdata",rootDir)

setList <- c("GBM")

colSets <- list(
	KIRC=c(rep("darkgreen",6),rep("orange",5),"red"),
	GBM=c(rep("darkgreen",5),rep("orange",4),"red")
)

# vertical separators
vLines <- list(
	GBM=c(5.5,9.5)
)

for (curSet in setList) {
	inFile <- sprintf("%s/%s_oneNetPer_FeatSel_results.Rdata",rootDir,curSet)
	load(inFile)
	
	colnames(val) <- sub("clinical","clin",colnames(val))
	mu <- colMeans(val,na.rm=TRUE)
	sem <- sapply(1:ncol(val),function(x) 
		sd(val[,x],na.rm=TRUE)/sqrt(nrow(val)))
	
	colSet <- colSets[[curSet]]
	plot(1:length(mu),mu,xlab="Data combinations",ylim=c(0.4,0.85),
		type='p',bty='n',ylab="AUCROC (mean+/-SEM)",xaxt='n',
		las=1,col=colSet,pch=16,cex=1.7)
	axis(1,at=1:length(mu), labels=colnames(val))
	segments(x0=1:length(mu), y0=mu-sem,y1=mu+sem,col=colSet,lwd=2)
	abline(h=0.5,col='red',lty=1,lwd=2)
	abline(h=0.7,col='red',lty=3,lwd=2)
	abline(v=vLines[[curSet]],lty=3,lwd=2)
	
	#boxplot(val,las=2,main=sprintf("%s (N=%i splits)",curSet,nrow(val)),
	#	bty='n',ylab="AUCROC")
	#abline(h=0.7,lty=3,col='red')
	#abline(h=0.5,lty=1,col='red')
	
	has_na <- colSums(is.na(val))
	if (any(has_na>0)) {
		cat("the following sections have 1+ failed rounds\n")
		print(has_na[which(has_na>0)])
	}
}
