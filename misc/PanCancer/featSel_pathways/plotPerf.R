#' Plot performance of various KIRC predictor conditions
#' Statistically compare them.

rm(list=ls())

# list of all conditions to collect data for and their i/o locations
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
outRoot <- sprintf("%s/2017_PanCancer_Survival",rootDir)

setFile <- "KIRCpathway_locations.txt"
setInfo	<- read.delim(setFile,sep="\t",h=T,as.is=T)
setInfo <- subset(setInfo, inc=="yes")

curSet <- "KIRC"

mega_roc <- matrix(NA,nrow=100,ncol=nrow(setInfo))
colnames(mega_roc) <- setInfo$name

mega_pr <- matrix(NA,nrow=100,ncol=nrow(setInfo))
colnames(mega_pr) <- setInfo$name

# treat clinical-one-net-only as its own feature
tmpFile <- sprintf("%s/%s", outRoot, setInfo$outdir[1])
lnames <- load(tmpFile)
mega_roc[,1] <- val[,1]
mega_pr[,1] <- val_pr[,1]
 
rm(val,val_pr)
for (ctr in 2:nrow(setInfo)) {
	cur <- setInfo$name[ctr]
	print(cur)
	maxk <- setInfo$maxK[ctr]
	outFile <- sprintf("%s/%s/KIRC_results_%s.Rdata",
		outRoot,setInfo$outdir[ctr],cur)
	
	lnames <- load(outFile)
	mega_roc[,ctr]	<- val; 
	mega_pr[,ctr]		<- val_pr
}

# same with single RNA net
tmpFile <- sprintf("%s/oneNetPer_FeatSel/KIRC_oneNetPer_FeatSel_results.Rdata", 
	outRoot)
load(tmpFile)
idx <- which(colnames(val)=="rna")
mega_roc <- cbind(mega_roc, val[,idx])
colnames(mega_roc)[ncol(mega_roc)] <- "rna"
rm(val)

tmpFile <- sprintf("%s/oneNetPer_FeatSel/AUCPR/KIRC_oneNetPer_FeatSel_results_prauc.Rdata",outRoot)
lnames<-load(tmpFile)
idx <- which(colnames(val)=="rna")
mega_pr <- cbind(mega_pr, val[,idx])
colnames(mega_pr)[ncol(mega_pr)] <- "rna"
rm(val)

#### one netPer
###inFile <- sprintf("%s/oneNetPer_FeatSel/%s_oneNetPer_FeatSel_results.Rdata",
###	rootDir,curSet)
###load(inFile)
###oneNetPer <- val; rm(val)
###
#### clinical by var & RNA by pathway
###inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_results.Rdata",
###	rootDir,curSet)
###load(inFile)
###oneNetPer <- cbind(oneNetPer,val[,1])
###colnames(oneNetPer)[ncol(oneNetPer)] <- "clinRNAPath"
###
#### clinical in one net & RNA by pathway
###inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_oneClinNet_results.Rdata",
###	rootDir,curSet)
###load(inFile)
###oneNetPer <- cbind(oneNetPer,val[,1])
###colnames(oneNetPer)[ncol(oneNetPer)] <- "clinOneRNAPath"
###
#### clinical one & RNA by pathway - RANDOMLY-sampled
###inFile <- sprintf("%s/randomNets_170508/KIRC_randomMean.Rdata",rootDir)
###load(inFile)
###oneNetPer <- cbind(oneNetPer, randomMean)
###colnames(oneNetPer)[ncol(oneNetPer)] <- "rnd"
###rm(randomMean)
###
####  RNA by pathway - RANDOMLY-sampled
###inFile <- sprintf("%s/pathOnly_randomNets_170503/KIRC_randomMean.Rdata",
###	rootDir)
###load(inFile)
###tmp <- c(randomMean, rep(NA, 100-length(randomMean)))
###oneNetPer <- cbind(oneNetPer, tmp)
###colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnlyrnd"
###print(length(randomMean))
###rm(randomMean)
###
#### clinical one & RNA by pathway - consensus nets
###inFile <- sprintf("%s/consensus_170502/KIRC_consensusRes.Rdata",rootDir)
###load(inFile)
###oneNetPer <- cbind(oneNetPer, consRes) 
###colnames(oneNetPer)[ncol(oneNetPer)] <- "consNet"
###
#### pathway only
###inFile <- sprintf("%s/pathwaysOnly_170502/KIRC_pathway_results.Rdata",rootDir)
###load(inFile)
###oneNetPer <- cbind(oneNetPer, consRes) 
###colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnly"
###
#### clinical by var only
###inFile <- sprintf("%s/clinNets_170430/%s_clinNets_170430_results.Rdata",
###	rootDir,curSet)
###load(inFile)
###oneNetPer <- cbind(oneNetPer,val[,1])
###colnames(oneNetPer)[ncol(oneNetPer)] <- "clinNets"
###
#### Pathway only - consensus nets
###inFile <- sprintf("%s/pathOnly_consensus_170509/KIRC_pathOnly_consensusRes.Rdata",
###	rootDir)
###load(inFile)
###oneNetPer <- cbind(oneNetPer, consRes) 
###colnames(oneNetPer)[ncol(oneNetPer)] <- "pathOnlycons"

colnames(mega_roc) <- sub("clinical","clin",colnames(mega_roc))

mega_roc <- mega_roc[,c(1:4,7,5:6)]
mega_pr <- mega_pr[,c(1:4,7,5:6)]

colSet <- c("red","red","purple","hotpink1", 
	"dodgerblue3","dodgerblue3","gray40")

#colSet <- c(rep("darkgreen",6),rep("orange",5),"red",
#		"purple","pink","blue","brown","gray","blue","pink")
#mega_roc <- mega_roc[,-c(2,4:11)]
#colSet		<- colSet[-c(2,4:11)]

###colSet <- c("red","purple","grey40","hotpink1",
###			"red","purple",
###			"dodgerblue3","dodgerblue3","grey40","hotpink1",
###			"green")

# do this at the end, after organizing mega_roc.

	postscript(sprintf("%s/KIRC_perf.eps",outRoot),width=18,height=6)
	tryCatch({
		par(bty='n',mar=c(3,5,2,4),mfrow=c(2,2))
for (cur_dat in c("roc","pr")) {
	print(cur_dat)
   if (cur_dat == "roc") curdat <- mega_roc
	 else curdat <- mega_pr

	mu <- colMeans(curdat,na.rm=TRUE)
	sem <- sapply(1:ncol(curdat),function(x) 
		sd(curdat[,x],na.rm=TRUE)/sqrt(nrow(curdat)))
	
	if (cur_dat =="roc") ylim <- c(0.65,0.9) else ylim <- c(0.6,0.85)
		
		lbl <- colnames(curdat)
		lbl[which(lbl=="clinPath")] <- "cOne\nRpath"
		lbl[which(lbl=="clinNetsPathNets")] <- "clin-nets\nRpath"
		lbl[which(lbl=="clinNets")] <- "clin-nets"
		# pathway effect
		idxSet <- list(pathways=5:7, clinPerf=1:4)
		for (nm in names(idxSet)) {
			idx <- idxSet[[nm]]
			plot(1:length(idx), mu[idx],ylim=ylim,
				type='n',bty='n',ylab="AUCROC\n(mean+/-SEM)",xaxt='n',
				las=1,cex.axis=1.4,xlim=c(1,length(idx)),cex.axis=1.4)
				abline(h=c(0.7,0.8),col='cadetblue3',lty=3,lwd=3)
				points(1:length(idx),mu[idx],type='p',col=colSet[idx],pch=16) 
					#,cex=0.5)
			title(curSet)
	
		# x-axis labels
			axis(1,at=1:length(idx), labels=lbl[idx],cex.axis=1)
	
		# error bars
		segments(x0=1:length(idx), y0=mu[idx]-sem[idx],
				y1=mu[idx]+sem[idx],col=colSet[idx],lwd=3)
		segments(x0=1:length(idx)-0.05, x1=1:length(idx)+0.05,
				y0=mu[idx]-sem[idx],y1=mu[idx]-sem[idx],col=colSet[idx],lwd=4)
		segments(x0=1:length(idx)-0.05, x1=1:length(idx)+0.05,
				y0=mu[idx]+sem[idx],y1=mu[idx]+sem[idx],col=colSet[idx],lwd=4)
		abline(h=0.5,col='red',lty=1,lwd=2)
	
		abline(v=c(6.5,12.5),col='black',lty=1)
	}	
	has_na <- colSums(is.na(curdat))
	if (any(has_na>0)) {
		cat("the following sections have 1+ failed rounds\n")
		print(has_na[which(has_na>0)])
	}
	
	
	# x,y = columns to compare
	# type = less,greater,two.sided. type of wmw test to run
	.wmwtest <- function(x,y,type) {
		wmw <- wilcox.test(curdat[,x],curdat[,y],alternative=type)
		cat(sprintf("%s vs %s\t\t%1.2f vs %1.2f\t\tWMW p < %1.2e\n",
			x,y, median(curdat[,x]), median(curdat[,y]),
			wmw$p.value))
	}
#if (cur_dat == "roc") {
	.wmwtest("clinOne","clinNets","less")
	.wmwtest("clinNets","clinNetsPathBest","less")
	.wmwtest("pathOnly","pathOnlyRnd","greater")
	.wmwtest("rna","pathOnly","less")
	###.wmwtest("pathOnlyrnd","pathOnly","less")
#	}
}
	}, error=function(ex){
		print(ex)
	}, finally={
		dev.off()
	})

	
