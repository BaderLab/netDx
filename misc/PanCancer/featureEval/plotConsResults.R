#' plot results of running predictor with only best consensus nets vs random

require(ROCR)

inDir <- "/mnt/data2/BaderLab/PanCancer_KIRC"

dirSet <- list(
		allCons=sprintf("%s/output/consNets_170410",inDir),
		randNets=sprintf("%s/output/randomNets_170410",inDir),
		bestCons=sprintf("%s/output/bestConsNets_170410",inDir)
	)

predSet <- list()
for (d in names(dirSet)) {
	fSet <- dir(sprintf("%s/predictions", dirSet[[d]]),pattern="prediction")
	cat(sprintf("%s: Got %i predictions\n", d, length(fSet)))
	val <- rep(NA, length(fSet))
	for (fName in fSet) {	
		dat <- read.delim(sprintf("%s/predictions/%s",dirSet[[d]],fName),
			sep="\t",h=T,as.is=T)
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
	predSet[[d]] <- val
}

