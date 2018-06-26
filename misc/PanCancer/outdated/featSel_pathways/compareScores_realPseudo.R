# compare net scores for real and pseudo pathways
# do real pathways score higher on average than pseudo?
rm(list=ls())

inDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/realPseudo_170810"

pctPass <- 0.7
for (gp in c("SURVIVEYES","SURVIVENO")) {
	cat("*********\n")
	cat(sprintf("%s\n",gp))
	cat("*********\n")
	
	inFile <- sprintf("%s/featSelNets/realPseudo_thresh10_pctPass0.70_%s_netScores.txt",inDir,gp)
	dat <- read.delim(inFile,sep="\t",h=T,as.is=T)
	idx <- grep("PSEUDO_", dat[,1])

	sets <- list(
		real=dat[setdiff(1:nrow(dat),idx),],
		pseudo=dat[idx,]	
	)

	out <- sapply(names(sets),function(nm) {
		cur <- sets[[nm]]
		cat(sprintf("%s\n", nm))
		maxNetS <- matrix(NA, nrow=length(cur[,1]),ncol=1)
		for (sc in 3:10) {
			tmp <- rowSums(cur[,-1] >= sc)
			idx <- which(tmp >= floor(pctPass * (ncol(cur)-1)))
			cat(sprintf("\t%i : %i pass\n", sc, length(idx)))
			maxNetS[idx,1] <- sc
	}
	idx <- which(!is.na(maxNetS))
	maxNetS <- maxNetS[idx,,drop=F]
	netNames <- cur[idx,1]

	#return(list(names=netNames,maxscore=maxNetS))
	cbind(nm, maxNetS)
})

	out <- as.data.frame(do.call("rbind",out))
	out[,1] <- as.character(out[,1])
	out[,2] <- as.numeric(as.character(out[,2]))

	require(ggplot2)
	p <- ggplot(out,aes(x=V2)) + geom_bar(aes(fill=as.factor(nm)),position=position_dodge())
	p <- p + ggtitle(sprintf("# nets passing score cutoffs >=70%% time\n%s",gp))
	print(p)
	browser()

}
