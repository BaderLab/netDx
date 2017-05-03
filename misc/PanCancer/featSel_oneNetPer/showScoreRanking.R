#' ranks features based on consistency across train/test splits
rm(list=ls())

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneNetPer_FeatSel/featSelNets"

require(plotrix)

tset <- 1:10
pset <- seq(0.25,1,0.25)

pdf(sprintf("%s/netScoreRankings.pdf",inDir),width=11,height=5)
tryCatch({
	par(mar=c(3,5,1,1),mfrow=c(2,1))
for (curSet in c("KIRC")) { #,"OV","LUSC","GBM")) {
	for (gp in c("SURVIVEYES","SURVIVENO")) {
		inFile <- sprintf("%s/%s_thresh10_%s_netScores.txt",inDir,curSet,gp)
		dat <- read.delim(inFile,sep="\t",h=T,as.is=T)
		netSet <- matrix(0,nrow=nrow(dat),ncol=10*4)
		rownames(netSet) <- dat[,1]
		
		ctr <- 1
		for (thresh in tset) {
			tmp <- rowSums(dat[,-1] >= thresh)
			for (pctPass in pset) {
				idx <- which(tmp >= floor(pctPass*ncol(dat[,-1])))
				netSet[idx,ctr] <- 1
				ctr <- ctr+1
			}
		}
	
		netSet <- netSet[order(rowSums(netSet),decreasing=TRUE),]
		
		plotrix::color2D.matplot(netSet,show.values=F,axes=F)
		axis(2,at=seq_len(nrow(dat))-0.5,
			labels=rev(sub(".profile|_cont","",
				rownames(netSet))),cex.axis=1,las=1,
			ylab="Feature",xlab="Cutoff and % pass") 
		axis(1,at=seq(1,40,4)-0.5,labels=tset,cex.axis=0.8)
		title(sprintf("%s:%s",curSet,gp))
	}
}},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

