#' ranks features based on consistency across train/test splits
rm(list=ls())

#inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/featSel_pathways_170426/featSelNets"
inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneClinNet_featSelNets"

require(plotrix)

tset <- 9:10
pset <- seq(0.1,1,0.1)

pctPass4file <- 1.0

pdf(sprintf("%s/netScoreRankings.pdf",inDir),width=11,height=5)
tryCatch({
	par(mar=c(2,30,1,1),mfrow=c(2,1))
for (curSet in c("KIRC")) { # ,"OV","LUSC","GBM")) {
	for (gp in c("SURVIVEYES","SURVIVENO")) {
		inFile <- sprintf("%s/%s__thresh10_pctPass%1.2f_%s_netScores.txt",
			inDir,curSet,pctPass4file,gp)
		dat <- read.delim(inFile,sep="\t",h=T,as.is=T)
		netSet <- matrix(0,nrow=nrow(dat),ncol=length(tset)*length(pset))
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
		netSet <- netSet[which(rowSums(netSet)>=1),]
		netSet <- 1-netSet
		
		plotrix::color2D.matplot(netSet,show.values=F,axes=F,xlab="",ylab="",
			extremes=c("blue","white"),border='gray50')
		axis(2,at=seq_len(nrow(netSet))-0.5,
			labels=rev(sub(".profile|_cont","",
				rownames(netSet))),cex.axis=0.6,las=1,
			ylab="Feature",xlab="Cutoff and % pass") 
		axis(1,at=seq(1,length(tset)*length(pset),length(pset))-0.5,
			labels=tset,cex.axis=0.8)
		title(sprintf("%s:%s",curSet,gp))
	}
}},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

