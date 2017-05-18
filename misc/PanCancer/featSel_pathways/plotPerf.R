#' Plot performance of various KIRC predictor conditions
#' Statistically compare them.

rm(list=ls())

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival"
outDir <- rootDir

setList <- c("KIRC")
for (curSet in setList) {
# one netPer
inFile <- sprintf("%s/oneNetPer_FeatSel/%s_oneNetPer_FeatSel_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- val; rm(val)

# clinical by var & RNA by pathway
inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinRNAPath"

# clinical in one net & RNA by pathway
inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_oneClinNet_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinOneRNAPath"

# clinical one & RNA by pathway - RANDOMLY-sampled
inFile <- sprintf("%s/randomNets_170508/KIRC_randomMean.Rdata",rootDir)
load(inFile)
oneNetPer <- cbind(oneNetPer, randomMean)
colnames(oneNetPer)[ncol(oneNetPer)] <- "rnd"
rm(randomMean)

#  RNA by pathway - RANDOMLY-sampled
inFile <- sprintf("%s/pathOnly_randomNets_170503/KIRC_randomMean.Rdata",
	rootDir)
load(inFile)
tmp <- c(randomMean, rep(NA, 100-length(randomMean)))
oneNetPer <- cbind(oneNetPer, tmp)
colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnlyrnd"
print(length(randomMean))
rm(randomMean)

# clinical one & RNA by pathway - consensus nets
inFile <- sprintf("%s/consensus_170502/KIRC_consensusRes.Rdata",rootDir)
load(inFile)
oneNetPer <- cbind(oneNetPer, consRes) 
colnames(oneNetPer)[ncol(oneNetPer)] <- "consNet"

# pathway only
inFile <- sprintf("%s/pathwaysOnly_170502/KIRC_pathway_results.Rdata",rootDir)
load(inFile)
oneNetPer <- cbind(oneNetPer, consRes) 
colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnly"

# clinical by var only
inFile <- sprintf("%s/clinNets_170430/%s_clinNets_170430_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinNets"

# Pathway only - consensus nets
inFile <- sprintf("%s/pathOnly_consensus_170509/KIRC_pathOnly_consensusRes.Rdata",
	rootDir)
load(inFile)
oneNetPer <- cbind(oneNetPer, consRes) 
colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnlycons"

colnames(oneNetPer) <- sub("clinical","clin",colnames(oneNetPer))

colSet <- c(rep("darkgreen",6),rep("orange",5),"red",
		"purple","pink","blue","brown","gray","blue","pink")
oneNetPer <- oneNetPer[,-c(2,4:11)]
colSet		<- colSet[-c(2,4:11)]

oneNetPer <- oneNetPer[,c("clin","clinOneRNAPath","rnd","consNet",
			"clinNets","clinRNAPath",
			"rna","pathOnly","pathOnlyrnd","pathOnlycons",
			"all")]
colSet <- c("red","purple","grey40","hotpink1",
			"red","purple",
			"dodgerblue3","dodgerblue3","grey40","hotpink1",
			"green")


# do this at the end, after organizing oneNetPer.
mu <- colMeans(oneNetPer,na.rm=TRUE)
sem <- sapply(1:ncol(oneNetPer),function(x) 
	sd(oneNetPer[,x],na.rm=TRUE)/sqrt(nrow(oneNetPer)))

postscript(sprintf("%s/KIRC_perf.eps",outDir),width=18,height=3)
tryCatch({
	par(bty='n',mar=c(3,5,2,4))

	plot(1:length(mu),mu,xlab="Data combinations",ylim=c(0.65,0.85),
		type='n',bty='n',ylab="AUCROC\n(mean+/-SEM)",xaxt='n',
		las=1,cex.axis=1.4)
	abline(h=c(0.7,0.8),col='cadetblue3',lty=3,lwd=3)
	abline(v=c(4.5,6.5,10.5),col='black')
	points(1:length(mu),mu,type='p',col=colSet,pch=16) #,cex=0.5)
	title(curSet)

	# x-axis labels
	lbl <- colnames(oneNetPer)
	lbl[which(lbl=="clinRNAPath")] <- "cNets\nRpath"
	lbl[which(lbl=="clinOneRNAPath")] <- "cOne\nRpath"
	lbl[which(lbl=="clinNets")] <- "clin-nets"
	axis(1,at=1:length(mu), labels=lbl,cex.axis=1.4)

	# error bars
	segments(x0=1:length(mu), y0=mu-sem,y1=mu+sem,col=colSet,lwd=3)
	segments(x0=(1:length(mu))-0.08, x1=(1:length(mu))+0.08,
			y0=mu-sem,y1=mu-sem,col=colSet,lwd=4)
	segments(x0=(1:length(mu))-0.08, x1=(1:length(mu))+0.08,
			y0=mu+sem,y1=mu+sem,col=colSet,lwd=4)
	abline(h=0.5,col='red',lty=1,lwd=2)

	#boxplot(oneNetPer,las=1,
	#	col=c(rep("antiquewhite",6), rep("chocolate1",6),
	#		"red","blue","purple"),
	#	main=sprintf("%s (N=%i splits)",curSet,nrow(oneNetPer)),
	#	bty='n',ylab="AUCROC",cex.axis=1.5,xaxt="n")

	abline(v=c(6.5,12.5),col='black',lty=1)
	
	#text(3.5,1,"Vanilla, single",cex=1.5)
	#text(9.5,1,"Vanilla, Clin+ 1 -omic",cex=1.5)
	#text(14.5,1, "Pathways/Vars",cex=1.5)
}, error=function(ex){
	print(ex)
}, finally={
	dev.off()
})

has_na <- colSums(is.na(oneNetPer))
if (any(has_na>0)) {
	cat("the following sections have 1+ failed rounds\n")
	print(has_na[which(has_na>0)])
}

###cat("-----\n")
###wmw <- wilcox.test(oneNetPer[,"clinArna"],oneNetPer[,"clinRNAPath"],
###	alternative="less")
###cat(sprintf("clinRNA vs clinRNAPath: %1.2f vs %1.2f; WMW (one-sided) p < %1.2e\n",
###	median(oneNetPer[,"clinArna"]), median(oneNetPer[,"clinRNAPath"]),
###	wmw$p.value))
###cat("-----\n")
###wmw <- wilcox.test(oneNetPer[,"clinArna"],oneNetPer[,"clinOneRNAPath"])
###cat(sprintf("clinRNA vs clinRNAPath: %1.2f vs %1.2f; WMW p < %1.2e\n",
###	median(oneNetPer[,"clinArna"]), median(oneNetPer[,"clinOneRNAPath"]),
###	wmw$p.value))
###cat("-----\n")
###wmw <- wilcox.test(oneNetPer[,"clin"],oneNetPer[,"clinArna"],
###			alternative="less")
###cat(sprintf("clin vs clinARNA: %1.2f vs %1.2f; WMW p < %1.2e\n",
###	median(oneNetPer[,"clin"]), median(oneNetPer[,"clinArna"]),
###	wmw$p.value))
###cat("-----\n")
###wmw <- wilcox.test(oneNetPer[,"clinRNAPath"],oneNetPer[,"clinNets"],
###			alternative="less")
###cat(sprintf("clinRNAPath vs clinNets: %1.2f vs %1.2f; WMW p < %1.2e\n",
###	median(oneNetPer[,"clinRNAPath"]), median(oneNetPer[,"clinNets"]),
###	wmw$p.value))
###cat("-----\n")

# x,y = columns to compare
# type = less,greater,two.sided. type of wmw test to run
.wmwtest <- function(x,y,type) {
	wmw <- wilcox.test(oneNetPer[,x],oneNetPer[,y],alternative=type)
	cat(sprintf("%s vs %s\t\t%1.2f vs %1.2f\t\tWMW p < %1.2e\n",
		x,y, median(oneNetPer[,x]), median(oneNetPer[,y]),
		wmw$p.value))
}

.wmwtest("clin","clinNets","less")
.wmwtest("clinOneRNAPath","rnd","greater")
.wmwtest("rna","pathOnly","less")
.wmwtest("pathOnlyrnd","pathOnly","less")

}

