# compare PRANK-based classification for real v shuffled data.
rm(list=ls())

require(ROCR)
source("rankPatients_test.R")

inDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/simRank_170719/realPathways"
shufDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/simRank_170719/designD_rmFSgenes_170717"
survFile <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/input/KIRC_binary_survival.txt"

pheno <- read.delim(survFile,sep="\t",h=T,as.is=T)
colnames(pheno)[1] <- c("ID")
survStr <- rep("SURVIVEYES",nrow(pheno))
survStr[which(pheno$is_alive<1)] <- "SURVIVENO"
pheno$STATUS <- survStr

# pFile1=yes, pFile2=no
.getOut <- function(pFile1,pFile2,pheno) {
real_yes <- rankPatient(pFile1)
real_no <- rankPatient(pFile2)
x <- merge(real_yes,real_no,by="ID")
x$GM_score.x <- order(x$GM_weight.x)/nrow(x)
x$GM_score.y <- order(x$GM_weight.y)/nrow(x)
x$GM_rank_yes <- rank(x$GM_weight.x)
x$GM_rank_no <- rank(x$GM_weight.y)
x$rankDiff <- x$GM_rank_yes - x$GM_rank_no

#x$SURVIVEYES_SCORE <- x$GM_rank_yes/nrow(x)
#x$SURVIVENO_SCORE <- x$GM_rank_no/nrow(x)
colnames(x)[which(colnames(x)=="GM_score.x")] <- "SURVIVEYES_SCORE"
colnames(x)[which(colnames(x)=="GM_score.y")] <- "SURVIVENO_SCORE"

predLbl <- rep(NA, nrow(x))
predLbl[x$rankDiff<0] <- "SURVIVEYES"
predLbl[x$rankDiff>1] <- "SURVIVENO"
x$PRED_LABEL <- predLbl

x <- merge(x,pheno,by="ID")
x
}

getConf <- function(dat,c1="SURVIVEYES") {
	dat <- na.omit(dat)
    tp <- sum(dat$STATUS==dat$PRED_LABEL & dat$STATUS == c1)
    tn <- sum(dat$STATUS==dat$PRED_LABEL & dat$STATUS != c1)
   fp <- sum(dat$STATUS!=dat$PRED_LABEL & dat$STATUS != c1)
    fn <- sum(dat$STATUS!=dat$PRED_LABEL & dat$STATUS == c1)

	return(c(tp,tn,fp,fn)/nrow(dat))
}

# ------------------------------------
# run work

source("getAUC.R")
realdat <- list()
shufdat <- list()
confmat <- matrix(ncol=8,nrow=50)
for (k in 2:50) {
	yesSfx <- sprintf("rng%i/SURVIVEYES/SURVIVEYES_query-results.report.txt.PRANK",k)
	noSfx <- sprintf("rng%i/SURVIVENO/SURVIVENO_query-results.report.txt.PRANK",k)
	
	pFile1 <- sprintf("%s/%s",inDir,yesSfx)
	pFile2 <- sprintf("%s/%s",inDir,noSfx)
	realres <- .getOut(pFile1,pFile2,pheno)
	tmp <- getAUC(realres)
	realdat[[k]] <- tmp
	confmat[k,1:4] <- getConf(realres)
	
	pFile1 <- sprintf("%s/%s",inDir,yesSfx)
	pFile1 <- sprintf("%s/%s",shufDir,yesSfx)
	pFile2 <- sprintf("%s/%s",shufDir,noSfx)
	shufres <- .getOut(pFile1,pFile2,pheno)
	tmp <- getAUC(shufres)
	shufdat[[k]] <- tmp
	confmat[k,5:8] <- getConf(shufres)
}

colnames(confmat) <- c(paste("real",c("tp","tn","fp","fn"),sep=""),
		paste("shuf",c("tp","tn","fp","fn"),sep=""))
confmat <- confmat[-1,]

require(reshape2)
df <- melt(confmat)
src <- rep("shuf",nrow(df))
src[grep("real",df$X2)] <- "real"
df$src <- factor(src,levels=c("real","shuf"))
df$X2 <- sub("real|shuf","",as.character(df$X2))
require(ggplot2)
p <- ggplot(df, aes(x=X2,y=value)) + geom_boxplot(aes(colour=src))

###par(mfrow=c(1,2))
###plot(0,0,xlim=c(0,1),ylim=c(0,1))
###for (k in 2:50) {
###	x <- realdat[[k]][["ROC"]]
###	lines(x[,1],x[,2],col=rgb(1,0,0,0.2))
###}
###
###for (k in 2:50) {
###	x <- shufdat[[k]][["ROC"]]
###	lines(x[,1],x[,2],col=rgb(0,0,0,0.2))
###}
###
