#' compare consensus nets for different scenarios
rm(list=ls())
require(ggplot2)
require(scales)
require(reshape2)

dSet <- list(
	pathOnly="/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/pathwaysOnly_170502/featSelNets",
	pathNoFS="/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/pathOnly_noFS_170809/featSelNets",
	pseudo_featSel="/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/pseudo_featSel_170809/featSelNets"
)
randomDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival/randomD_pseudoPath_noPathGenes_170804"

pctPass <- 0.7

for (gp in c("SURVIVEYES","SURVIVENO")) {
	gpNets <- sapply(names(dSet),function(nm) {
		x <- dSet[[nm]]
		cat(sprintf("%s\n",x))
		fName <- dir(x,pattern=sprintf("%s_netScores.txt", gp))
		netS <- read.delim(sprintf("%s/%s",x,fName),sep="\t",h=T,as.is=T)
		netNames <- netS[,1]; netS <- netS[,-1]
	
	# compute the max score per net for pctPass % of trials
	maxNetS <- matrix(NA, nrow=length(netNames),ncol=1)
	for (sc in 1:10) {
			tmp <- rowSums(netS >= sc)
			idx <- which(tmp >= floor(pctPass * ncol(netS)))
			cat(sprintf("\t%i : %i pass\n", sc, length(idx)))
			maxNetS[idx,1] <- sc
	}
	idx <- which(is.na(maxNetS[,1]))
	if (any(idx)) maxNetS[idx,1] <- 0
	cbind(nm,netNames,maxNetS)
	})

	df <- data.frame(do.call("rbind",gpNets))
	df[,1] <- as.character(df[,1])
	df[,2] <- as.character(df[,2])
	df[,3] <- as.numeric(as.character(df[,3]))
	colnames(df) <- c("setting","netNames","best_score")
	#df <- subset(df, best_score>=3)
	
	require(ggplot2)
	p <- ggplot(df,aes(best_score)) + geom_density(aes(fill=setting))
	p <- p + ggtitle(gp)
	
	p <- ggplot(df,aes(best_score)) + geom_bar(aes(fill=setting),position="dodge")# + scale_y_continuous(labels=percent_format()
	p <- p+ggtitle(gp)
	#print(p)

}	

# compare % of splits for which a net scores 9 or 10 in real vs
# randomly-sampled data

par(las=1,bty='n',mfrow=c(1,2),mar=c(2,4,4,2))
for (gp in c("SURVIVEYES","SURVIVENO")) {
		nm <- "pseudo_featSel"
		x <- dSet[[nm]]
		cat(sprintf("%s\n",x))
		fName <- dir(x,pattern=sprintf("%s_netScores.txt", gp))
		netS <- read.delim(sprintf("%s/%s",x,fName),sep="\t",h=T,as.is=T)
		netNames <- netS[,1]; netS <- netS[,-1]
		netPass <- rowSums(netS >=9,na.rm=T)
		realdat <- data.frame(netName=sub(".profile","",netNames),
			real_netScore=netPass)
		realdat[,1] <- as.character(realdat[,1])
		
		# read random data
		rnd <- read.delim(sprintf("%s/%s.netTally.txt",randomDir,gp),sep="\t",
			h=T,as.is=T)
		colnames(rnd)[2] <- c("random_netScore")
	x <- merge(x=realdat,y=rnd,by="netName")
	y <- list(real=x[,2],pseudo_pathways=x[,3])
	names(y)[1] <- nm
	boxplot(y,ylab="# times a net scores >=9 (out of 100)",
		main=sprintf("%s\n # (out of 100) net scores >=9",gp))
}
