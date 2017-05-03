#' collect results for oneNetPer.
require(ROCR)

dirList <- list(
	KIRC="/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_oneNetPer_170426",
	LUSC="/home/netdx/BaderLab/PanCancer_LUSC/output/featSel_oneNetPer_170425", # VM4
	OV="/home/ahmad/tcga_datasets/OV/output/featSel_oneNetPer_170425",# VM2
	GBM="/home/spai/BaderLab/PanCancer_GBM/output/featSel_oneNetPer_170425" #VM3
)
outDirList <- list(
	KIRC="/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival",
	LUSC="/home/netdx/BaderLab/PanCancer_common", 	# VM4
	OV="/home/spai/BaderLab/PanCancer_common", # VM2
	GBM="/home/spai/BaderLab/PanCancer_common" # VM3
)

# iterations that have run so far
maxK <- list(
	KIRC=100,
	LUSC=21,
	GBM=19,
	OV=40
)

outList <- list()
megaList <- list()
combList <- list(    
    clinical="clinical.profile",    
	mir="mir.profile",
	rna="rna.profile",
	prot="prot.profile",
	cnv="cnv.profile",
	dnam="dnam.profile",
    clinicalArna=c("clinical.profile","rna.profile"),    
    clinicalAmir=c("clinical.profile","mir.profile"),    
    clinicalAprot=c("clinical.profile","rppa.profile"),    
    clinicalAdnam=c("clinical.profile","dnam.profile"),    
    clinicalAcnv=c("clinical.profile","cnv.profile"),    
    all="all")  

for (cur in "KIRC") { #names(dirList)) {
	curd <- dirList[[cur]]
	currCombList <- combList
	if (cur=="LUSC") {# no dnam
		currCombList[["dnam"]] <- NULL
		currCombList[["clinicalAdnam"]] <- NULL
	} else if (cur == "GBM") {
		currCombList[["prot"]] <- NULL
		currCombList[["clinicalAprot"]] <- NULL
	}
	maxk <- maxK[[cur]]
	kset <- 1:maxk
	cat(sprintf("Num runs=%i\n", maxk))

		val <- matrix(NA,nrow=length(kset),ncol=length(currCombList))
		colnames(val) <- names(currCombList)
		for (k in kset) {
			for (nm in names(currCombList)) {
				finDir	<-sprintf("%s/rng%i/%s",curd,k,nm)
				dat <- read.delim(sprintf("%s/predictionResults.txt",finDir),
						,sep="\t",h=T,as.is=T)
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
				idx <- which(colnames(val)==nm)
				val[k,idx] <- performance(pred, "auc")@y.values[[1]]
				}
			}
		}

	outFile <- sprintf("%s/%s_oneNetPer_FeatSel_results.Rdata",
		outDirList[[cur]],cur)
	cat("Saving to file.\n")
	save(val,file=outFile)
}
