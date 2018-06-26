#' ranks features based on consistency across train/test splits
rm(list=ls())

inDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab/2017_PanCancer_Survival/clinRNA_best"

require(plotrix)

tset <- 10
pset <- seq(0.1,1,0.1)

pctPass4file <- 0.7

pdf(sprintf("%s/netScoreRankings.pdf",inDir),width=8,height=5)
tryCatch({
	par(mar=c(2,33,1,1),mfrow=c(2,1))
for (curSet in c("KIRC")) { # ,"OV","LUSC","GBM")) {
	for (gp in c("SURVIVEYES","SURVIVENO")) {
		inFile <- sprintf("%s/clinRNA_best_thresh10_pctPass%1.2f_%s_netScores.txt",
			inDir,pctPass4file,gp)
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

		tmp <- tolower(rev(sub(".profile|_cont","",
				rownames(netSet))))
		tmp <- gsub("_"," ",tmp)	

		axis(2,at=seq_len(nrow(netSet))-0.5,
			labels=tmp,cex.axis=0.6,las=1,
			ylab="Feature",xlab="Cutoff and % pass") 
		axis(1,at=seq(1,length(tset)*length(pset),length(pset))-0.5,
			labels=tset,cex.axis=0.8)
		title(sprintf("%s:%s",curSet,gp),cex.main=0.8)
	}
}},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

