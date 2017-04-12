#' plot results of running predictor with only best consensus nets vs random

require(ROCR)

rootDir <- "/mnt/data2/BaderLab" # VM1
#rootDir <- "/home/netdx/BaderLab" # VM4

dt <- format(Sys.Date(),"%y%m%d")
saveData <- TRUE # set to true to save, false to plot

if (saveData) {
for (curSet in c("KIRC")) {
	inDir <- sprintf("%s/PanCancer_%s",rootDir,curSet)
	dirs <- dir(path=sprintf("%s/output", inDir),pattern="randomNets")
	dirSet <- list()
	for (d in dirs) dirSet[[d]] <- sprintf("%s/output/%s",inDir,d)
	
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
		}
		predSet[[d]] <- val
	}
	pdf(sprintf("%s/PanCancer_common/%s_consRes_plot_random.pdf",
		rootDir,curSet),width=8,height=3)
	tryCatch({
		layout(matrix(c(1,1,2),ncol=3,nrow=1,byrow=TRUE))
		boxplot(predSet,main=sprintf("%s: %i iterations", 
			curSet,length(predSet[[1]])),ylab="AUC")
		tmp <-unlist(lapply(predSet,mean))
		boxplot(tmp,
			main="Mean over 25 resamplings\n(10 random samples)",
			ylab="mean AUCROC over 25 resamplings",ylim=c(0.4,1))
		abline(h=0.5,lty=3,col='red')
		cat("Summary of 10 random resamplings\n")
		print(summary(tmp))
	},error=function(ex){
		print(ex)
	},finally={
		dev.off()
	})

	#save(predSet, file=sprintf("%s/PanCancer_common/%s_consRes_%s.Rdata",rootDir,curSet,dt))
}
} else { # load and plot
	inDir <- sprintf("%s/PanCancer_common", rootDir)
	res <- list()
	for (curSet in c("KIRC","LUSC","OV","GBM")) {
		fName <- dir(path=inDir,pattern="Rdata")
		fName <- fName[grep(curSet, fName)]
		load(sprintf("%s/%s",inDir,fName))
		res[[curSet]] <- do.call("rbind",predSet)
	}

	blah <- list()
	for (nm in names(res)) { 
		x <- res[[nm]]
		df <- as.data.frame(x)
		df$type <- rownames(df)
		df$cancertype <- nm
		blah[[nm]] <- as.matrix(df)
	}
	blah2 <- do.call("rbind",blah)
	blah3 <- suppressWarnings(data.frame(blah2, stringsAsFactors=FALSE))
	for (k in 1:25) blah3[,k] <- as.numeric(blah3[,k])
	
	require(reshape2)
	blah4 <- melt(blah3)
	
	require(ggplot2)                                                             
	p <- ggplot(blah4, aes(cancertype, value)) + geom_boxplot(aes(colour=type))  
	p <- p + ylab("AUCROC (25 runs)")
	p <- p + ggtitle("PanCancer test perf on consensus/random nets") 
	p <- p + geom_hline(aes(yintercept=0.5),colour='red')
	                                                                             
	pdfFile <- sprintf("%s/consRes_perf_%s.pdf", inDir,dt)                       
	pdf(pdfFile,width=8,height=4)                                                
	print(p)                                                                     
	dev.off() 
}
