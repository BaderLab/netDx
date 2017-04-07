# plot AUCROC of all four cancer types
rm(list=ls()) 
require(ggplot2)
require(ROCR)

dirBase <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	GBM=sprintf("%s/2017_TCGA_GBM/output/ownTrain_170206",dirBase),
	OV=sprintf("%s/2017_TCGA_OV/output/ownTrain_170205",dirBase),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/ownTrain_170205",dirBase),
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/netdx_integrate_170206",dirBase)
	)

dt <- format(Sys.Date(),"%y%m%d")
featResFile <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/featSel_perf_170227.Rdata"

sem <- function(x) sd(x)/sqrt(length(x))

outList <- list()
megaList <- list()
full2 <- list()
full <- list()
for (cur in names(dirList)) {
	print(cur)

	subdList <- dir(path=sprintf("%s/run1",dirList[[cur]]))
	for (subd in subdList) { # per datatype
		cat(sprintf("\t%s\n",subd))
		val <- rep(NA,100)
		for (k in 1:100) {
				dat <- read.delim(sprintf("%s/run%i/%s/predictionResults.txt",
							dirList[[cur]],k,subd),sep="\t",h=T,as.is=T)
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
				curdat <- data.frame(cancer=cur,
																 datatype=subd,
																 AUC=val,
												 				stringsAsFactors=FALSE) 
		avg <- aggregate(curdat$AUC, 
					by=list(cancer=curdat$cancer, datatype=curdat$datatype),
					FUN=mean);
		colnames(avg)[3] <- "mean"
		myse <- aggregate(curdat$AUC, 
					by=list(cancer=curdat$cancer, datatype=curdat$datatype),
					FUN=sem)
		colnames(myse)[3] <- "sem"
		x <- merge(x=avg,y=myse,by=c("cancer","datatype"))
		outList[[cur]][[subd]]  <- x

		# unit-level data
		full2[[cur]][[subd]] <- curdat
	}
	megaList[[cur]]	<- do.call("rbind",outList[[cur]])
	full[[cur]]			<- do.call("rbind", full2[[cur]])
}

out <- do.call("rbind",megaList)

# add feature selection as its own 'datatype'
lnames <- load(featResFile)
featSel_agg$datatype <- "pathway_features"
out <- rbind(out, featSel_agg)
featSel_full <- do.call("rbind",featSel_full)
browser()
featSel_full <- cbind(featSel_full,datatype="pathway_features")
featSel_full <- featSel_full[,c(1,3,2)]

out$datatype <- factor(out$datatype, 
										 levels=c("clinical","clinicalArna","clinicalAmir",
															"clinicalArppa","all","pathway_features"))

colList <- list(
			clinical=rgb(178,24,43,max=255),
			clinicalArna=rgb(249,116,36,max=255),
			clinicalAmir=rgb(95,168,204,max=255),
			clinicalArppa=rgb(237,66,32,max=255),
			all="grey50",
			pathway_features="purple")

# add eb
limits <- aes(ymax=mean+sem,ymin=mean-sem)

# plot and save pdf
p <- ggplot(out,aes(colour=datatype,y=mean,x=cancer))
p <- p + geom_point(position=position_dodge(width=0.9),size=2.5) +
			geom_errorbar(position=position_dodge(width=0.9),limits,width=0.4)
p <- p + scale_colour_manual(values=unlist(colList))
p <- p + ylim(0.5,0.85)
# b/w theme, larger axis tick labels
p <- p + theme_bw() + theme(axis.text=element_text(size=14)) 
		
p <- p + ggtitle("netDx: PanCancer survival AUCROC")
pdf(sprintf("%s/PanCancer_featSel_survival_%s.pdf",dirname(featResFile),
						dt),width=10,height=3)
print(p);
dev.off()

# compare LUSC clinical to clinical+ RPPA
alldat <- do.call("rbind",full)
alldat <- rbind(alldat, featSel_full)

idx1 <- which(alldat$cancer %in% "LUSC" & alldat$datatype %in% "clinical")
idx2 <- which(alldat$cancer %in% "LUSC" & alldat$datatype %in% "clinicalArppa")
wmw<- wilcox.test(alldat[idx1,"AUC"],alldat[idx2,"AUC"],
									alternative="less")
cat(sprintf("LUSC: clinical = %1.2f ; clin + RPPA = %1.2f (WMW p < %1.2e)\n",
						mean(alldat[idx1,"AUC"]),mean(alldat[idx2,"AUC"]),wmw$p.value))

# compare KIRC clinical to clinical+ RNA
idx1 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "clinical")
idx2 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "clinicalArna")
wmw<- wilcox.test(alldat[idx1,"AUC"],alldat[idx2,"AUC"],
		alternative="less")
cat(sprintf("KIRC: clinical = %1.2f ; clin + RNA = %1.2f (WMW p < %1.2e)\n",
						mean(alldat[idx1,"AUC"]),mean(alldat[idx2,"AUC"]),wmw$p.value))

# compare KIRC best to KIRC pathways
idx1 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "clinicalArna")
idx2 <- which(alldat$cancer %in% "KIRC" & alldat$datatype %in% "pathway_features")
wmw<- wilcox.test(alldat[idx1,"AUC"],alldat[idx2,"AUC"],
		alternative="less")
cat(sprintf("KIRC: clinArna = %1.2f ; pathway = %1.2f %1.2f (WMW p < %1.2e)\n",
						mean(alldat[idx1,"AUC"]),mean(alldat[idx2,"AUC"]),
						mean(alldat[idx1,"AUC"])-mean(alldat[idx2,"AUC"]),
						wmw$p.value))

