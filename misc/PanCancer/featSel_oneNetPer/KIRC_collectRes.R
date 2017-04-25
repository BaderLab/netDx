#' collect results for oneNetPer.

dirList <- list(
	KIRC="/mnt/data2/BaderLab/PanCancer_KIRC/output/featSel_oneNetPer_170424"
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

for (cur in names(dirList)) {
	curd <- dirList[[cur]]
	maxk <- 1
	kset <- 1:maxk
	cat(sprintf("Num runs=%i\n", maxk))

		val <- matrix(NA,nrow=length(kset),ncol=length(combList))
		colnames(val) <- names(combList)
		for (k in kset) {
			for (nm in c("clinicalArna","prot","all")) {
				finDir	<-sprintf("%s/rng%i/%s",curd,k,nm)
				dat <- read.delim(sprintf("%s/predictionResults.txt",finDir),
						,sep="\t",h=T,as.is=T)
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
