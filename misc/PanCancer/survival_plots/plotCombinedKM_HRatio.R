# average KM curves for KIRC 
rm(list=ls())

require(rms)
require(survival)
require(survminer)
source("plot.survfit.custom.R")

tumourType <- "KIRC"

rootDir <-  "/Users/shraddhapai/Dropbox/netDx/BaderLab"
if (tumourType=="KIRC") {
	dataDir <- sprintf("%s/2017_TCGA_KIRC/output/pruned_180204",rootDir)
	survFile <- sprintf("%s/2017_TCGA_KIRC/input/KIRC_OS_core.txt",rootDir)
	clinFile <- sprintf("%s/2017_TCGA_KIRC/input/KIRC_clinical_core.txt",rootDir)
} else if (tumourType=="LUSC") {
	dataDir <- sprintf("%s/2017_TCGA_LUSC/output/LUSC_oneNetPer_170425",rootDir)
	survFile <- sprintf("%s/2017_TCGA_LUSC/input/LUSC_OS_core.txt",rootDir)
	clinFile <- sprintf("%s/2017_TCGA_LUSC/input/LUSC_clinical_core.txt",rootDir)
} else if (tumourType=="OV") {
	dataDir <- sprintf("%s/2017_TCGA_OV/output/OV_oneNetPer_170425",rootDir)
	survFile <- sprintf("%s/2017_TCGA_OV/input/OV_OS_core.txt",rootDir)
	clinFile <- sprintf("%s/2017_TCGA_OV/input/OV_clinical_core.txt",rootDir)
}

survDat <- read.delim(survFile,sep="\t",h=T,as.is=T)
clinDat <- read.delim(clinFile,sep="\t",h=T,as.is=T)
pheno <- merge(x=survDat,y=clinDat,by="feature")

plotDF <- list() # compiles survival curves across all iterations
megaDF <- list()
hratio <- c() # cum hazards ratio for all iterations
for (k in 1:100) {
	print(k)
	dat <- read.delim(sprintf("%s/rng%i/clinical/predictionResults.txt",
		dataDir,k),sep="\t",h=T,as.is=T)
	colnames(dat)[1] <- "feature"
	dat <- merge(x=dat,y=pheno,by="feature")

	# force first entry to be YES and second to be NO so we can tell them apart
	# in the output.
	dat$PRED_CLASS <- factor(dat$PRED_CLASS,levels=c("SURVIVEYES","SURVIVENO"))

	megaDF[[k]] <- dat

 	dat$SurvObj <- with(dat, Surv(OS_OS, STATUS_INT == 0))

	# get cum hazards ratio for this split
	model <- coxph(SurvObj~PRED_CLASS, data=dat)
	hratio <- c(hratio,summary(model)$coef[1,2])

	fit <- npsurv(SurvObj ~ PRED_CLASS, data = dat,
		conf.type = "log-log")

	par(mfrow=c(1,2))
	out <- plot.survfit.custom(fit)
	plot(0,0,type='n',xlim=c(0,max(out$ends$x)),ylim=c(0,1))
	out[[1]] <- as.data.frame(out[[1]])
	out[[2]] <- as.data.frame(out[[2]])
#	lines(out[[1]]$xx,out[[1]]$yy,col='green')
#	lines(out[[2]]$xx,out[[2]]$yy,col='red')

	newdf <- out[[1]]; newdf$PRED_CLASS <- "SURVIVEYES"; newdf$split <- k
	newdf2 <- out[[2]]; newdf2$PRED_CLASS <- "SURVIVENO"; newdf2$split <- k
	
	plotDF[[k]] <- rbind(newdf,newdf2)
}

# approach 1: pool all results and make a single KM curve
res <- do.call("rbind",megaDF)
res$SurvObj <- with(res, Surv(OS_OS,STATUS_INT==0))
fit <- npsurv(SurvObj ~ PRED_CLASS, data=res,conf.type="log-log")

idx <- which(hratio > 50)
if (any(idx)) hratio <- hratio[-idx]
hratio <- data.frame(group="tumour",hratio=hratio)
p <- ggplot(hratio,aes(group,y=hratio)) + geom_boxplot() + ylim(c(0,quantile(hratio$hratio,0.98))+3)
p <- p + ggtitle(sprintf("Cum hazard ratio (one per split)(N=%i)",nrow(hratio)))
p <- p + geom_hline(yintercept=1,lty=2)
p <- p + theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

p1 <- survminer::ggsurvplot(fit,data=res,conf.int=TRUE)
p2 <- p

# approach 2 : pool all step functions
res <- do.call("rbind",plotDF)
# compute ci
out <- list()
for (k in unique(res$PRED_CLASS)) {
	res2 <- subset(res,PRED_CLASS==k)
	ub <- c(); lb <- c(); xx <- unique(res2$xx); muy <- c()
	for (x in xx) {
		yy <- res2$yy[which(res2$xx == x)]
		mu <- mean(yy); offset <- sd(yy)/sqrt(length(yy))
		lb <- c(lb, mu-offset)
		ub <- c(ub, mu+offset)
		muy <- c(muy, mu)
	}
	out[[k]] <- data.frame(x=xx,y=muy,lb=lb,ub=ub,PRED_CLASS=k)
}
blah <- do.call("rbind",out)
p3 <- ggplot(blah,(aes(x=x,y=y,colour=PRED_CLASS))) + geom_line() + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.2) + ggtitle("Manual mean+CI of compiled KM") + ylab("% survival") +xlab("time (months)")

pdf(sprintf("%s.pdf",tumourType))
print(p1); print(p2); print(p3)
dev.off()
