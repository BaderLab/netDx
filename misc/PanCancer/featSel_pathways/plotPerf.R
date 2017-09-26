#' Plot performance of various KIRC predictor conditions
#' Statistically compare them.
rm(list=ls())

# list of all conditions to collect data for and their i/o locations
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
outRoot <- sprintf("%s/2017_PanCancer_Survival",rootDir)
dt <- format(Sys.Date(),"%y%m%d")

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
	if (cur == "rnaOne") {
		outFile <- sprintf("%s/oneNetPer_FeatSel/fromAhmad_170428/all_rdata/KIRC_oneNetPer_FeatSel_results.Rdata",outRoot)
		lnames <- load(outFile)
		mega_roc[,ctr] <- val[,which(colnames(val)=="rna")]; rm(val)
		outFile <- sprintf("%s/oneNetPer_FeatSel/AUCPR/KIRC_oneNetPer_FeatSel_results_prauc.Rdata",outRoot)
		lnames <- load(outFile)
		mega_pr[,ctr] <- val[,which(colnames(val)=="rna")]; rm(val)
	} else {
		maxk <- setInfo$maxK[ctr]
		if (cur %in% "pathOnlyRnd_old") {
		outFile <- sprintf("%s/%s/KIRC_results_pathOnlyRnd.Rdata",
			outRoot,setInfo$outdir[ctr],cur)
		} else if (cur %in% c("pathOnlyRnd_noFS","pathOnlyRnd_Shuf","pathOnlyRnd_25")) {
		outFile <- sprintf("%s/%s/KIRC_randomMean.Rdata",
			outRoot,setInfo$outdir[ctr])
		} else {
		outFile <- sprintf("%s/%s/KIRC_results_%s.Rdata",
			outRoot,setInfo$outdir[ctr],cur)
		}
	
		lnames <- load(outFile)
		mega_roc[,ctr]	<- val; 
		mega_pr[,ctr]		<- val_pr
	} 
}

colnames(mega_roc) <- sub("clinical","clin",colnames(mega_roc))
keepnm <- c("clinOne","clinNets","clinNetsPathBest",
		"rnaOne","pathOnly","pseudo_featSel",
	"scrambled2")
	#	"pathFull_AltClass", "pathAlt_noFS",
	#	"pathOnly80","pathOnly90","pathOnly95")
mega_roc <- mega_roc[,which(colnames(mega_roc) %in% keepnm)]
mega_pr <- mega_pr[,which(colnames(mega_pr) %in% keepnm)]

postscript(sprintf("%s/KIRC_perf_%s.eps",outRoot,dt),width=6,height=6)
	tryCatch({
		par(bty='n',mar=c(3,4,1,1),mfrow=c(2,1))#mfrow=c(2,2))
for (cur_dat in c("roc","pr")) {
	print(cur_dat)
   if (cur_dat == "roc") curdat <- mega_roc
	 else curdat <- mega_pr

	mu <- colMeans(curdat,na.rm=TRUE)
	sem <- sapply(1:ncol(curdat),function(x) 
		sd(curdat[,x],na.rm=TRUE)/sqrt(nrow(curdat)))
	
	if (cur_dat =="roc") {
		ylim <- c(0.45,0.9) 
		ylab <- "AUROC(mean+/-SEM)"
	} else {
		ylim <- c(0.45,0.85)
		ylab <- "AUPR(mean+/-SEM)"
	}
		
		lbl <- colnames(curdat)
		lbl[which(lbl=="clinPath")] <- "cOne\nRpath"
		lbl[which(lbl=="clinNets")] <- "clin-nets"

		# pathway effect
		idxSet <- list(pathways=4:length(mu), clinPerf=1:3,
				AltClass=c(4,9:11))
		#	pathways=c("pathOnly","pathRnd_D_noFS",
		#		"pathOnly80","pathOnly90","pathOnly95",
		#		"pathFull_AltClass","pathAlt_noFS")
		#	clinPerf=c("clinOne","clinNets","clinNetsPathBest"))
		for (nm in names(idxSet)[1]) {
			idx <- idxSet[[nm]]
			plot(1:length(idx), mu[idx],ylim=ylim,
				type='n',bty='n',ylab=ylab,xaxt='n',cex.lab=1,xlab="",
				las=1,cex.axis=0.7,xlim=c(0.5,length(idx)+0.5),
				srt=45)

				abline(h=c(0.7,0.8),col='cadetblue3',lty=3,lwd=3)
				clrs <- numeric()
				for (i in names(mu)[idx]) {
					clrs <- c(clrs, setInfo$color[which(setInfo$name==i)])
				}
				points(1:length(idx),mu[idx],type='p',col=clrs,pch=16) 

			ln <- colSums(!is.na(curdat[,idx]))
			text(x=1:length(idx), y=0.6,labels=sprintf("N=%i",ln))

	
		# x-axis labels
			axis(1,at=1:length(idx), labels=sub("pathOnlyRnd_","",lbl[idx]),
				cex.axis=0.7)
	
		# error bars
		segments(x0=1:length(idx), y0=mu[idx]-sem[idx],
				y1=mu[idx]+sem[idx],col=clrs,lwd=3)
		segments(x0=1:length(idx)-0.05, x1=1:length(idx)+0.05,
				y0=mu[idx]-sem[idx],y1=mu[idx]-sem[idx],col=clrs,lwd=4)
		segments(x0=1:length(idx)-0.05, x1=1:length(idx)+0.05,
				y0=mu[idx]+sem[idx],y1=mu[idx]+sem[idx],col=clrs,lwd=4)
		abline(h=0.5,col='red',lty=1,lwd=2)
	
		#abline(v=c(6.5,12.5),col='black',lty=1)
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
			x,y, median(curdat[,x],na.rm=TRUE), median(curdat[,y],na.rm=TRUE),
			wmw$p.value))
	}
#if (cur_dat == "roc") {
	cat(sprintf("Stat tests for %s\n------------------\n",cur_dat))
#	.wmwtest("clinOne","clinNets","less")
#	.wmwtest("clinNets","clinNetsPathBest","less")
#	.wmwtest("pathOnly","pathOnlyRnd","greater")
#	.wmwtest("pathOnly","pathOnlyCons","less")
#	.wmwtest("pathOnly","pathRnd_D_noFS","greater")
	#.wmwtest("pathOnly","pathFull_AltClass","less")
#	.wmwtest("pathOnly","pathOnly80","less")
#	.wmwtest("pathOnly","pathOnly90","less")
#	.wmwtest("pathOnly","pathOnly95","less")
	#.wmwtest("pathFull_AltClass","pathAlt_noFS","greater")
	.wmwtest("rnaOne","pathOnly","less")
#	.wmwtest("pathOnly","pathOnly_noFS","greater")
#	}
}
	}, error=function(ex){
		print(ex)
	}, finally={
	dev.off()
	})
