#' collect results for clinNetsOnly
require(ROCR)

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/KIRC_clinNets_170430",rootDir),
	 # VM
	LUSC="",
	OV="",
	GBM=""
)
outDirList <- list(
	KIRC=sprintf("%s/2017_PanCancer_survival",rootDir), 	# VM1
	LUSC="",#"/home/netdx/BaderLab/PanCancer_common", 	# VM4
	OV="",#"/home/spai/BaderLab/PanCancer_common", # VM2
	GBM=""#"/home/spai/BaderLab/PanCancer_common" # VM3
)

# iterations that have run so far
maxK <- list(
	KIRC=100
	#LUSC=21,
	#GBM=19,
	#OV=40
)

outList <- list()
megaList <- list()

for (cur in "KIRC") { 
	curd <- dirList[[cur]]

	maxk <- maxK[[cur]]
	kset <- 1:maxk
	cat(sprintf("Num runs=%i\n", maxk))

		val <- matrix(NA,nrow=maxk,ncol=1)
		colnames(val) <- ""
		for (k in kset) {
		print(k)
					finDir	<-sprintf("%s/rng%i",curd,k)
				dat <- read.delim(sprintf("%s/predictionResults.txt",finDir),
						sep="\t",h=T,as.is=T)
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
					val[k,1] <- performance(pred, "auc")@y.values[[1]]
			}
		}

	outFile <- sprintf("%s/%s_clinNets_170430_results.Rdata",
		outDirList[[cur]],cur)
	cat("Saving to file.\n")
	save(val,file=outFile)
}
