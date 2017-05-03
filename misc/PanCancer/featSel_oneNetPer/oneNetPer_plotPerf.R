# plot AUCROC for predictor
rm(list=ls())
require(reshape2)
require(multcomp)

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneNetPer_FeatSel"

inDir <- sprintf("%s/fromAhmad_170428/all_rdata",rootDir)

setList <- c("KIRC","GBM","LUSC","OV")

colSets <- list(
	KIRC=c(rep("darkgreen",6),rep("orange",5),"red"),
	GBM=c(rep("darkgreen",5),rep("orange",4),"red"),
	LUSC=c(rep("darkgreen",5),rep("orange",4),"red"),
	OV=c(rep("darkgreen",6),rep("orange",5),"red")
)

# vertical separators
vLines <- list(
	GBM=c(5.5,9.5),
	KIRC=c(6.5,11.5),
	LUSC=c(5.5,9.5),
	OV=c(6.5,11.5)
)

logFile <- sprintf("%s/perf.log.txt",rootDir)

pdf(sprintf("%s/perf.pdf",rootDir),width=13,height=6)
sink(logFile,split=TRUE)
tryCatch({
par(mfrow=c(2,2),mar=c(2,4,1,0))
for (curSet in setList) {
	cat("---------------------\n")
	cat(sprintf("%s\n",curSet))
	cat("---------------------\n")
	inFile <- sprintf("%s/%s_oneNetPer_FeatSel_results.Rdata",rootDir,curSet)
	load(inFile)
	
	colnames(val) <- sub("clinical","clin",colnames(val))
	mu <- colMeans(val,na.rm=TRUE)
	sem <- sapply(1:ncol(val),function(x) 
		sd(val[,x],na.rm=TRUE)/sqrt(nrow(val)))
	
	colSet <- colSets[[curSet]]
	plot(1:length(mu),mu,xlab="Data combinations",ylim=c(0.4,0.85),
		type='p',bty='n',ylab="AUCROC (mean+/-SEM)",xaxt='n',
		las=1,col=colSet,pch=16,cex=1.5,cex.axis=1.4)
	title(curSet)

	axis(1,at=1:length(mu), labels=colnames(val),cex.axis=1.4)
	segments(x0=1:length(mu), y0=mu-sem,y1=mu+sem,col=colSet,lwd=3)
	abline(h=0.5,col='red',lty=1,lwd=2)
	abline(h=0.7,col='red',lty=3,lwd=2)
	abline(v=vLines[[curSet]],lty=1,lwd=2)
	
	#boxplot(val,las=2,main=sprintf("%s (N=%i splits)",curSet,nrow(val)),
	#	bty='n',ylab="AUCROC")
	#abline(h=0.7,lty=3,col='red')
	#abline(h=0.5,lty=1,col='red')lbl	
	has_na <- colSums(is.na(val))
	if (any(has_na>0)) {
		cat("the following sections have 1+ failed rounds\n")
		print(has_na[which(has_na>0)])
	}
	
	# run ANOVA for effect of integration
	x <- melt(val)
	group <- rep("single",nrow(x))
	group[which(x$Var2 %in% c("clinArna","clinAmir","clinAdnam","clinAprot",
			"clinAcnv","all"))] <- "integrated"
	x$group <- as.factor(group)

	fit <- summary(aov(value~group,data=x))
#print(fit)
	cat(sprintf("P(F-test), integration changes effect = %1.2e\n", 
		fit[[1]][["Pr(>F)"]][1]))

	x_full <- x; rm(x)

set.seed(103)
	cat("\t------\n")
	# integration > clinical
	cat("B. Dunnett's test: clinical > integrated:\n")
	cat("\n")
	x <- subset(x_full, Var2 %in% c("clin", "clinArna","clinAmir","clinAdnam",
			"clinAprot","clinAcnv","all"))
	x$Var2 <- as.character(droplevels(x$Var2))
	x$Var2[which(x$Var2=="clin")] <- "A_clinical" # Dunnett compares
	x$Var2 <- as.factor(x$Var2)

			# first level to others
	aov_res <- aov(value ~ Var2, data = x)
	dunn <- summary(glht(aov_res, linfct=mcp(Var2="Dunnett")))
	t_neg <- sum(dunn$test$tstats<0) == length(dunn$test$tstats)
	psig	<- sum(dunn$test$pvalues < 0.05) == length(dunn$test$pvalues)
	if (t_neg & psig) {
		cat("\tClinical BETTER than integrated\n")
	} else {
 		cat("\tClinical NOT better than integrated\n")
	}	
	print(dunn)
	cat("\t------\n")

	if (curSet %in% c("GBM","OV")) {
	cat("C. Dunnett's test: clinical > other single\n")
cat("\n")
	x <- subset(x_full, Var2 %in% c("clin", "rna","mir","dnam",
			"prot","cnv"))
	x$Var2 <- as.character(droplevels(x$Var2))
	x$Var2[which(x$Var2=="clin")] <- "A_clinical" # Dunnett compares
	x$Var2 <- as.factor(x$Var2)

			# first level to others
	aov_res <- aov(value ~ Var2, data = x)
	dunn <- summary(glht(aov_res, linfct=mcp(Var2="Dunnett")))
	t_neg <- sum(dunn$test$tstats<0) == length(dunn$test$tstats)
	psig	<- sum(dunn$test$pvalues < 0.05) == length(dunn$test$pvalues)
	if (t_neg & psig) {
		cat(sprintf("%s: Clinical BETTER than other single\n",curSet))
	} else {
 		cat(sprintf("%s: Clinical NOT better than other single\n",curSet))
	}	
	print(dunn)
	} else if (curSet=="LUSC") {
		cat("C. Dunnett's test: proteomics > other single\n")
		cat("\n")
		x <- subset(x_full, Var2 %in% c("clin", "rna","mir","dnam",
				"prot","cnv"))
		x$Var2 <- as.character(droplevels(x$Var2))
		# Dunnett compares first to others
		x$Var2[which(x$Var2=="prot")] <- "A_prot" 
		x$Var2 <- as.factor(x$Var2)
	
				# first level to others
		aov_res <- aov(value ~ Var2, data = x)
		dunn <- summary(glht(aov_res, linfct=mcp(Var2="Dunnett")))
		t_neg <- sum(dunn$test$tstats<0) == length(dunn$test$tstats)
		psig	<- sum(dunn$test$pvalues < 0.05) == length(dunn$test$pvalues)
		if (t_neg & psig) {
			cat(sprintf("%s: Proteomics BETTER than other single\n",curSet))
		} else {
	 		cat(sprintf("%s: Protemics NOT better than other single\n",curSet))
		}	
		print(dunn)

		cat(sprintf("%s: WMW: Proteomic != clinical\n", curSet))
		cat(sprintf("p < %1.2e\n",
			wilcox.test(val[,"clin"], val[,"prot"])$p.value))
	} else if (curSet=="KIRC") {
		cat("C. Dunnett's test: clin+DNAm > other paired\n")
		cat("\n")
		x <- subset(x_full, Var2 %in% c("clinAdnam", "clinArna","clinAmir",
				"clinAprot","clinAcnv"))
		x$Var2 <- as.character(droplevels(x$Var2))
		# Dunnett compares first to others
		x$Var2[which(x$Var2=="clinAdnam")] <- "A_clinDNAm" 
		x$Var2 <- as.factor(x$Var2)
	
				# first level to others
		aov_res <- aov(value ~ Var2, data = x)
		dunn <- summary(glht(aov_res, linfct=mcp(Var2="Dunnett")))
		t_neg <- sum(dunn$test$tstats<0) == length(dunn$test$tstats)
		psig	<- sum(dunn$test$pvalues < 0.05) == length(dunn$test$pvalues)
		if (t_neg & psig) {
			cat(sprintf("%s: Clin+DNAm BETTER than other integrated\n",curSet))
		} else {
	 		cat(sprintf("%s: Clin+DNAm NOT better than other integrated\n",curSet))
		}	
		print(dunn)
	}

}
},error=function(ex){
	print(ex)
},finally= {
	dev.off()
	sink(NULL)
})
