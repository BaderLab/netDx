# compare ranking of test patients in real vs shuffled pathways
rm(list=ls())

realDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/pathway_170502"
shufDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/pathRandom_designD_rmFSgenes_170717"

survFile <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/input/KIRC_binary_survival.txt"

outDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_PanCancer_Survival"
dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/simRank_%s",outDir,dt)

pheno <- read.delim(survFile,sep="\t",h=T,as.is=T)

source("compareTestSimRanks.R")
realDat <- compareTestSimRanks(realDir)
x <- as.numeric(realDat)
x <- x[!is.na(x)]

shufDat <- compareTestSimRanks(shufDir)
y <- as.numeric(shufDat)
y <- y[!is.na(y)]

midx <- match(rownames(realDat),rownames(shufDat))
if (all.equal(rownames(shufDat)[midx],rownames(realDat))!=TRUE) {
	cat("dont match"); browser()
}
shufDat <- shufDat[midx,]

pdf(sprintf("%s/group_YESNOdiff_violins.pdf",outDir))
require(caroline)
violins(list(real=x,shuf=y),ylab="(YES-NO)_score")
title("Delta similarity of test patients")
dev.off()

df	<- cbind("real",x)
df2 <- cbind("shuf",y)
df	<- as.data.frame(rbind(df,df2))
df[,2] <- as.numeric(as.character(df[,2]))
df[,3] <- abs(df[,2])
colnames(df) <- c("type","simBias","abs_simBias")

require(ggplot2)
source("multiplot.R")

ks <- ks.test(df[which(df[,1]=="real"),3], df[which(df[,1]=="shuf"),3])
p <- ggplot(df, aes(simBias,fill=type,colour=type)) 
p <- p + geom_density(alpha=0.1)
p <- p + ggtitle(sprintf("simBias (YES_rank-NO_rank)\n"))
p <- p + geom_vline(xintercept=0)

p2 <- ggplot(df, aes(abs_simBias,fill=type,colour=type)) + geom_density(alpha=0.1) 
p2 <- p2 + ggtitle(sprintf("abs(YES_rank-NO_rank)\nKS p < %1.2e",ks$p.value)) 
p2 <- p2 + geom_vline(xintercept=0)

pdf(sprintf("%s/group_YESNO_density.pdf",outDir),width=11,height=5)
multiplot(plotlist=list(p,p2),cols=2)
dev.off()

# ---------------
# now look at individual rankings
require(reshape2)
x <- na.omit(melt(t(realDat)))
x <- cbind(x,type="real")
y <- na.omit(melt(t(shufDat)))
y <- cbind(y,type="random")
z <- rbind(x,y)
z$type <- factor(as.character(z$type),levels=c("random","real"))
z$X2 <- as.character(z$X2)
z <- merge(x=z,y=pheno,by.x="X2",by.y="feature")

# get a difference in similarity bias for each patient
df <- numeric()
for (k in unique(z$X2)) {
	cur <- subset(z,X2 == k)
	cur_real <- mean(cur$value[cur$type=="real"])
	cur_random <- mean(cur$value[cur$type=="random"])
	df <- c(df, cur_real-cur_random)
}

df <- data.frame(name=unique(z$X2),df=df)
df <- merge(x=df,y=pheno,by.x="name",by.y="feature")
df$is_alive <- factor(df$is_alive)
p <- ggplot(df,aes(x=is_alive,y=df))+geom_boxplot() 
p <- p + ylab("patient level difference in (YES-NO)\n(real - random)")
p <- p + ggtitle("(YES-NO)_real - (YES-NO)_random\nDifference in rankings")
p <- p + geom_hline(yintercept=0)
pdf(sprintf("%s/patientLevel_MeanRankingDiff.pdf",outDir))
print(p)
dev.off()

# compare real and fake assignments for a random sample of 10 patients
pid <- sample(rownames(realDat),20,F)
z2 <- subset(z, X2 %in% pid)
p <- ggplot(z2, aes(x=X2,y=value)) + geom_boxplot(aes(colour=type)) 
p <- p + geom_hline(yintercept=0) + ylab("YES_score-NO_score") 
p <- p + scale_color_manual(values=c("#aaaaaa","#cc0000"))
pdf(sprintf("%s/patientLevel_RankingDiff.pdf",outDir),width=11,height=6)
print(p)
dev.off()

