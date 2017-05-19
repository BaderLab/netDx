#' collect results for oneNetPer.
#' clinical affected only, with normDiff
require(ROCR)

dirList <- list(
	KIRC="/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_oneNetPer_normDiff_170518",
	GBM="/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/output/GBM_oneNetPer_normDiff_170518")

outDir="/Users/shraddhapai/Documents/Research/BaderLab/2017_PanCancer_Survival/oneNetPer_FeatSel"

# iterations that have run so far
maxK <- list(
	KIRC=70,
	GBM=96
)


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
# --------------------------------------------------------------
# work begins
outList <- list()
megaList <- list()
combList <- list(    
    clinical="clinical_cont",    
    clinicalArna=c("clinical_cont","rna.profile"),    
    clinicalAmir=c("clinical_cont","mir.profile"),    
    clinicalAprot=c("clinical_cont","rppa.profile"),    
    clinicalAdnam=c("clinical_cont","dnam.profile"),    
    clinicalAcnv=c("clinical_cont","cnv.profile"),    
    all="all")  

for (cur in "GBM") { #names(dirList)) {
	curd <- dirList[[cur]]
	currCombList <- combList
	if (cur == "GBM") {
		currCombList[["prot"]] <- NULL
		currCombList[["clinicalAprot"]] <- NULL
	}
	maxk <- maxK[[cur]]
	kset <- 1:maxk
	cat(sprintf("Num runs=%i\n", maxk))

		val <- matrix(NA,nrow=length(kset),ncol=length(currCombList))
		val_pr<- matrix(NA,nrow=length(kset),ncol=length(currCombList))
		colnames(val) <- names(currCombList)
		colnames(val_pr) <- names(currCombList)
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
				# precision-recall
				tmp <- performance(pred,"prec","rec")
				val_pr[k,idx] <- prauc(tmp)
				}
			}
		}

	outFile <- sprintf("%s/%s_oneNetPer_clin_normDiff.Rdata",
		outDir,cur)
	cat("Saving to file.\n")
	save(val,val_pr,file=outFile)
}
