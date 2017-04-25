
rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneNetPer_featSel"

setList <- c("GBM")

for (curSet in setList) {
inFile <- sprintf("%s/%s_oneNetPer_FeatSel_results.Rdata",rootDir,curSet)
load(inFile)

colnames(val) <- sub("clinical","clin",colnames(val))
boxplot(val,las=2,main=sprintf("%s (N=%i splits)",curSet,nrow(val)),
	bty='n',ylab="AUCROC")
abline(h=0.7,lty=3,col='red')
abline(h=0.5,lty=1,col='red')

has_na <- colSums(is.na(val))
if (any(has_na>0)) {
	cat("the following sections have 1+ failed rounds\n")
	print(has_na[which(has_na>0)])
}
}
