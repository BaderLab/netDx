# plot AUCROC of all four cancer types
rm(list=ls()) 
require(ggplot2)
require(ROCR)

dirBase <- "/Users/shraddhapai/Documents/Research/BaderLab"

# with feature selection
dirList <- list(
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/output/featSel_incMutRPPA_170223",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222",dirBase)
	)

sem <- function(x) sd(x,na.rm=T)/sqrt(sum(!is.na(x)))

outList <- list()
megaList <- list()
full2 <- list()
full <- list()
for (cur in names(dirList)) {
	print(cur)

	if (cur == "GBM") {
				maxk <- 78
				kset <- 1:maxk
	} else if (cur == "OV") {
				maxk <- 48
				kset <- setdiff(1:maxk,40)
	} else {
				maxk <- 100
				kset <- 1:maxk
	}
	cat(sprintf("Num runs=%i\n", maxk))

		val <- rep(NA,maxk)
		for (k in kset) {
				dat <- read.delim(sprintf("%s/rng%i/predictionResults.txt",
							dirList[[cur]],k),sep="\t",h=T,as.is=T)
				pred <- prediction(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE,
						  dat$STATUS=="SURVIVEYES")

				c1 <- "SURVIVEYES" #numc[1]
				tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
				tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
				fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
				fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)
		
				# dummy score column needed for perfCalc()
				tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn) 
			  val[k] <- performance(pred, "auc")@y.values[[1]]
			}
		curdat <- data.frame(cancer=cur,AUC=val,stringsAsFactors=FALSE) 
		avg <- aggregate(curdat$AUC, 
				 by=list(cancer=curdat$cancer),FUN=mean,na.rm=T);
		colnames(avg)[2] <- "mean"
		myse <- aggregate(curdat$AUC, 
				 by=list(cancer=curdat$cancer),FUN=sem)
		colnames(myse)[2] <- "sem"
		x <- merge(x=avg,y=myse,by=c("cancer"))

		# unit-level data
	megaList[[cur]]	<- x
	full[[cur]]	 <- curdat
}

out <- do.call("rbind",megaList)

cat("---------------------------\n")
cat("Average performance\n")
print(out)
cat("---------------------------\n")
cat("Best performance\n")
x <- sapply(names(full), function(nm) 
						cat(sprintf("%s: %1.2f\n", nm, max(full[[nm]][,2],na.rm=T))))

# add eb
limits <- aes(ymax=mean+sem,ymin=mean-sem)

# plot and save pdf
p <- ggplot(out,aes(y=mean,x=cancer))
p <- p + geom_point(position=position_dodge(width=0.9),size=2.5) +
			geom_errorbar(position=position_dodge(width=0.9),limits,width=0.4)
p <- p + scale_colour_manual(values=unlist(colList))
# b/w theme, larger axis tick labels
p <- p + theme_bw() + theme(axis.text=element_text(size=14)) 
browser()
		
p <- p + ggtitle("netDx: PanCancer survival AUCROC")
pdf("PanCancer_basic_survival.pdf",width=10,height=3)
print(p);
dev.off()

# compare LUSC clinical to clinical+ RPPA
alldat <- do.call("rbind",full)
idx1 <- which(alldat$cancer %in% "LUSC" & alldat$datatype %in% "clinical")
idx2 <- which(alldat$cancer %in% "LUSC" & alldat$datatype %in% "clinicalArppa")
wmw<- wilcox.test(alldat[idx1,"AUC"],alldat[idx2,"AUC"],
									alternative="less")
cat(sprintf("LUSC: clinical = %1.2f ; clin + RPPA = %1.2f (WMW p < %1.2e)\n",
						mean(alldat[idx1,"AUC"]),mean(alldat[idx2,"AUC"]),wmw$p.value))

# compare KIRC clinical to clinical+ RNA
alldat <- do.call("rbind",full)
idx1 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "clinical")
idx2 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "clinicalArna")
wmw<- wilcox.test(alldat[idx1,"AUC"],alldat[idx2,"AUC"],
		alternative="less")
cat(sprintf("KIRC: clinical = %1.2f ; clin + RNA = %1.2f (WMW p < %1.2e)\n",
						mean(alldat[idx1,"AUC"]),mean(alldat[idx2,"AUC"]),wmw$p.value))
