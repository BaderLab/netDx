# see the net score overlap

inDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/consNetsTest"

for (curSet in "GBM") {
	maxScore <- NULL
	for (cutoff in 6:9) {
		dat <- read.delim(sprintf("%s/%s_thresh%i_SURVIVEYES_consNets.txt",
			inDir,curSet,cutoff),sep="\t",h=T,as.is=T)
		if (is.null(maxScore)) {
			maxScore <- data.frame(NET=dat[,1],SCORE=cutoff)
		} else {
			idx <- which(dat[,1] %in% maxScore[,1])
			maxScore[idx,2] <- cutoff
			idx <- which(!dat[,1] %in% maxScore[,1])
			if (any(idx)) {
			cat(sprintf("** warn: found %i not in maxScore\n",length(idx)))
			browser()
			}
		}
	}
	maxScore <- maxScore[order(maxScore[,2],decreasing=TRUE),]
###	require(plotrix)
###pdf("test.pdf",width=11,height=8)
###	par(mar=c(1,40,1,1))
###	plotrix::color2D.matplot(as.matrix(maxScore[,2]),
###		show.values=TRUE,axes=F,
###		xlab="",ylab="",vcex=1)
###	axis(2,at=seq_len(nrow(maxScore))-0.5,
###			labels=sub(".profile","",
###				sub("_cont","",rev(maxScore[,1]))),tick=F,
###			las=1,cex.axis=0.8)
###	dev.off()
	colnames(maxScore) <- c("net","max_score")
	write.table(maxScore,file=sprintf("%s/%s_SURVIVEYES_consNet_maxScore.txt",
		inDir,curSet),sep="\t",col=T,row=F,quote=F)
}
