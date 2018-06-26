#' plot BRCA results

require(netDx)
require(netDx.examples)
data(TCGA_BRCA)

rootDir <- "/home/shraddhapai/BaderLab/2017_BRCA/output"

pheno$STATUS[which(pheno$STATUS!="LumA")] <- "other"

# get one rna net perf
inDir <- sprintf("%s/BRCA_OneRNAnet_180221/.",rootDir)
outDir <- sprintf("%s/BRCA_OneRNAnet_180221/plot",rootDir)
predClasses <- unique(pheno$STATUS)
postscript(sprintf("%s/perf.eps",outDir))
predPerf_oneRNA <- plotPerf(inDir, predClasses=predClasses)
dev.off()

inDir_path <- sprintf("%s/BRCA_part2_180223", rootDir)
outDir <- sprintf("%s/BRCA_part2_180223/plot",rootDir)
postscript(sprintf("%s/perf.eps",outDir))
predPerf_path <- plotPerf(inDir_path, predClasses=predClasses)
dev.off()

auc_onerna <- unlist(lapply(predPerf_oneRNA,function(x) x$auroc))
auc_path <- unlist(lapply(predPerf_path,function(x) x$auroc))
wmw <- wilcox.test(auc_onerna,auc_path)
pdf("BRCA_oneNetVsPath.pdf")
boxplot(list(OneNet=auc_onerna, Pathways=auc_path),
	main=sprintf("BRCA, one vs path\n(WMW p < %1.2e)",
		wmw$p.value))
dev.off()
cat(sprintf("One RNA: AUC=%1.2f +/- %1.2f\n",mean(auc_onerna),sd(auc_onerna)))
cat(sprintf("Pathways: AUC=%1.2f +/- %1.2f\n",mean(auc_path),sd(auc_path)))


