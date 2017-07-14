#' plot results of running predictor with only best consensus nets vs random
rm(list=ls())

require(ROCR)

rootDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab"
inDir <- sprintf("%s/2017_TCGA_KIRC/output/KIRC_pathOnly_randomShuf_170712",
		rootDir) # VM1
outDir	<- sprintf("%s/2017_PanCancer_Survival/pathOnly_randomShuf_170712",
		rootDir)
curSet	<- "KIRC"
outFile	<- sprintf("%s/%s_randomMean.Rdata",outDir,curSet) 

if (!file.exists(outDir)) dir.create(outDir)

dt <- format(Sys.Date(),"%y%m%d")
saveData <- TRUE # set to true to save, false to plot

if (saveData) {
	dirs <- dir(path=inDir,pattern="randomNets")
	dirSet <- list()
	for (d in dirs) dirSet[[d]] <- sprintf("%s/%s",inDir,d)
	
	predSet <- list()
	for (d in names(dirSet)) {
		fSet <- dir(sprintf("%s/predictions", dirSet[[d]]),
			pattern="prediction")
		cat(sprintf("%s: Got %i predictions\n", d, length(fSet)))
		val <- rep(NA, length(fSet))
		ctr <- 1
		for (fName in fSet) {	
			dat <- read.delim(sprintf("%s/predictions/%s",dirSet[[d]],fName),
				sep="\t",h=T,as.is=T)
			if (nrow(dat)>1) {
	    	pred <- prediction(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE,    
	                          dat$STATUS=="SURVIVEYES") 

			c1 <- "SURVIVEYES" #numc[1]                                     
			tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)        
			tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)        
			fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)        
			fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)        
	
			# dummy score column needed for perfCalc()                      
			tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn)              
			val[ctr] <- performance(pred, "auc")@y.values[[1]]        
			ctr <- ctr+1
			} else {
				cat(sprintf("\t%s:%s: too few records!\n",d,fName))
			}
		}
		predSet[[d]] <- na.omit(val)
	}

} else { # load and plot
	load(outFile)
} 
	pdf(sprintf("%s/%s_consRes_plot_random.pdf",
		outDir,curSet),width=8,height=3)
	tryCatch({
		layout(matrix(c(1,1,2),ncol=3,nrow=1,byrow=TRUE))
		boxplot(predSet,main=sprintf("%s: %i iterations", 
			curSet,length(predSet[[1]])),ylab="AUC")
		tmp <-unlist(lapply(predSet,mean))
	
		# distribution of mean AUROC across all random resamplings
		boxplot(tmp,
			main=sprintf("Mean over 25 resamplings\n(%i random samples)",
				length(predSet)),
			ylab="mean AUROC over 25 resamplings",ylim=c(0.4,1))
		abline(h=c(0.5,0.7),lty=3,col='red')
		cat(sprintf("Summary of %i random resamplings\n",length(tmp)))
		print(summary(tmp))

		randomMean <- tmp
		save(randomMean,file=outFile)

	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})



