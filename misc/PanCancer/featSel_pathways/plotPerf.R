
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival"
outDir <- rootDir

setList <- c("KIRC")

for (curSet in setList) {
# one netPer
inFile <- sprintf("%s/oneNetPer_FeatSel/%s_oneNetPer_FeatSel_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- val; rm(val)

# clinical by var & RNA by pathway
inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinRNAPath"

# clinical in one net & RNA by pathway
inFile <- sprintf("%s/featSel_pathways_170426/%s_pathway_oneClinNet_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinOneRNAPath"

# clinical by var only
inFile <- sprintf("%s/clinNets_170430/%s_clinNets_170430_results.Rdata",
	rootDir,curSet)
load(inFile)
oneNetPer <- cbind(oneNetPer,val[,1])
colnames(oneNetPer)[ncol(oneNetPer)] <- "clinNets"

colnames(oneNetPer) <- sub("clinical","clin",colnames(oneNetPer))
#pdf(sprintf("%s/%s_perf.pdf", outDir,curSet), width=11,height=5)
tryCatch({
	par(bty='n',mar=c(3,4,2,0))
	boxplot(oneNetPer,las=1,
		col=c(rep("antiquewhite",6), rep("chocolate1",6),
			"red","blue","purple"),
		main=sprintf("%s (N=%i splits)",curSet,nrow(oneNetPer)),
		bty='n',ylab="AUCROC",cex.axis=1.5,xaxt="n")
	lbl <- colnames(oneNetPer)
	lbl <- sub("A","\n",lbl)
	lbl[which(lbl=="clinRNAPath")] <- "clin-nets\nRNA-path"
	lbl[which(lbl=="clinOneRNAPath")] <- "clinOne\nRNA-path"
	lbl[which(lbl=="clinNets")] <- "clin-nets"
	axis(1,at=1:ncol(oneNetPer), labels=lbl,cex.axis=1)

	abline(v=c(6.5,12.5),col='black',lty=1)
	abline(h=0.7,lty=3,col='red')
	abline(h=0.5,lty=1,col='red')
	
	text(3.5,1,"Vanilla, single",cex=1.5)
	text(9.5,1,"Vanilla, Clin+ 1 -omic",cex=1.5)
	text(14.5,1, "Pathways/Vars",cex=1.5)
}, error=function(ex){
	print(ex)
}, finally={
#	dev.off()
})

has_na <- colSums(is.na(oneNetPer))
if (any(has_na>0)) {
	cat("the following sections have 1+ failed rounds\n")
	print(has_na[which(has_na>0)])
}

cat("-----\n")
wmw <- wilcox.test(oneNetPer[,"clinArna"],oneNetPer[,"clinRNAPath"],
	alternative="less")
cat(sprintf("clinRNA vs clinRNAPath: %1.2f vs %1.2f; WMW (one-sided) p < %1.2e\n",
	median(oneNetPer[,"clinArna"]), median(oneNetPer[,"clinRNAPath"]),
	wmw$p.value))
cat("-----\n")
wmw <- wilcox.test(oneNetPer[,"clinArna"],oneNetPer[,"clinOneRNAPath"])
cat(sprintf("clinRNA vs clinRNAPath: %1.2f vs %1.2f; WMW p < %1.2e\n",
	median(oneNetPer[,"clinArna"]), median(oneNetPer[,"clinOneRNAPath"]),
	wmw$p.value))
cat("-----\n")
wmw <- wilcox.test(oneNetPer[,"clin"],oneNetPer[,"clinArna"],
			alternative="less")
cat(sprintf("clin vs clinARNA: %1.2f vs %1.2f; WMW p < %1.2e\n",
	median(oneNetPer[,"clin"]), median(oneNetPer[,"clinArna"]),
	wmw$p.value))
cat("-----\n")
wmw <- wilcox.test(oneNetPer[,"clinRNAPath"],oneNetPer[,"clinNets"],
			alternative="less")
cat(sprintf("clinRNAPath vs clinNets: %1.2f vs %1.2f; WMW p < %1.2e\n",
	median(oneNetPer[,"clinRNAPath"]), median(oneNetPer[,"clinNets"]),
	wmw$p.value))

}

