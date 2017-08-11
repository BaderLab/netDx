#' collect AUCROC and AUCPR for all conditions pathway related
require(ROCR)

# --------------------------------------------------------------
# helper routines

# given output of performance("precall") compute AUC-PR
prauc <- function(dat) {
	x <- dat@x.values[[1]] # recall
	y <- dat@y.values[[1]] # precision

	# remove NAN
	idx <- which(is.nan(y))
	if (any(idx)) { x <- x[-idx]; y <- y[-idx]}

	#x <- c(0,x,1); y <- c(1,y,0) # make sure points go from 0,0 to 1,1

	pracma::trapz(x,y)
}

# gets AUCROC and AUCPR from a table
.getAUC <- function(curFile) {
dat <- read.delim(curFile,sep="\t",h=T,as.is=T)
aucroc <- NA; aucpr <- NA
if (nrow(dat)>0) {
	pred <- prediction(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE,
	dat$STATUS=="SURVIVEYES")
	
	c1 <- "SURVIVEYES" #numc[1]
	tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
	tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
	fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
	fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)
	
	# dummy score column needed for perfCalc()
	tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn) 
	aucroc <- performance(pred, "auc")@y.values[[1]]
	# precision-recall
	tmp <- performance(pred,"prec","rec")
	aucpr <- prauc(tmp)
}
return(c(aucroc,aucpr))
}

# --------------------------------------------------------------
# list of all conditions to collect data for and their i/o locations
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
inRoot <- sprintf("%s/2017_TCGA_KIRC/output",rootDir)
outRoot <- sprintf("%s/2017_PanCancer_Survival",rootDir)

# --------------------------------------------------------------
dt <- format(Sys.Date(),"%y%m%d")

outList <- list()
megaList <- list()

inDirRoot <- sprintf("%s/pathSize_170808",inRoot)

pDir <- dir(inDirRoot,"pSize")
out <- list()
maxk <- 25
for (curp in pDir) {
	gDir <- dir(sprintf("%s/%s", inDirRoot,curp),pattern="numG")
	gDir <- sprintf("%s/%s/%s",inDirRoot,curp,gDir)
	
	for (curg in gDir) {
		cat(sprintf("%s:%s:Num runs=%i\n", curp,basename(gDir),length(kset)))
	
		rngDirs <- dir(curg,pattern="rng")
		val		<- matrix(NA,nrow=length(rngDirs),ncol=1)
		val_pr <- matrix(NA,nrow=length(rngDirs),ncol=1)
ctr <- 1
		for (k in rngDirs) {
				print(k)
				curFile<-sprintf("%s/%s/predictionResults.txt",curg,k)
				tmp <- .getAUC(curFile);
				val[ctr,1] <- tmp[1]
				val_pr[ctr,1] <- tmp[2]
	ctr <- ctr+1
		}
		baseG <- basename(curg)
		out[[sprintf("%s_%s",curp,baseG)]] <- list(AUROC=val,AUPR=val_pr)
	}
}

# order by pathway size
nm <- names(out)
upos <- regexpr("_",nm)
nm <- sub("pSize","",substr(nm, 1,upos-1))
idx <- order(as.integer(nm))

out <- out[idx]

roc <- unlist(lapply(out,function(x) mean(x[[1]])))
names(roc) <- names(out)
roc_sem <- unlist(lapply(out,function(x) sd(x[[1]])/sqrt(length(x[[1]]))))
pr <- unlist(lapply(out,function(x) mean(x[[2]])))
pr_sem <- unlist(lapply(out,function(x) sd(x[[2]])/sqrt(length(x[[2]]))))

par(mfrow=c(2,1),mar=c(4,4,2,2),las=1,cex.axis=1.1,bty='n')
plot(1:length(roc),roc,ylim=c(0.5,0.8),xaxt='n',
		xlab="Num pathways sampled for GM db for test classification\nEach pathway contains 10 randomly-sampled non-pathway genes (pseudo pathway)",cex=1.2,pch=16)
segments(x0=1:length(roc), y0=roc-roc_sem,
	y1=roc+roc_sem,lwd=3)
axis(side=1,at=1:length(roc),labels=names(out))
title("ROC")
abline(h=seq(0.6,0.8,0.05),lty=3,col='grey50')
# mark N
ln <- unlist(lapply(out,function(x) length(x[[1]])))
text(1:length(roc),0.55, sprintf("N=%i",ln))

plot(1:length(pr),pr,ylim=c(0.5,0.8),xaxt='n',cex=1.2,pch=16)
axis(side=1,at=1:length(pr),labels=names(out))
segments(x0=1:length(pr), y0=pr-pr_sem,
	y1=pr+pr_sem,lwd=3)
title("PR")
abline(h=seq(0.6,0.8,0.05),lty=3,col='grey50')
