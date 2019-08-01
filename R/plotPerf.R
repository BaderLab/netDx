#' Plots various measures of predictor performance for binary classifiers
#'
#' @details Plots individual and average ROC/PR curves. mean+/-SEM performance for a predictor run using nested
#' cross-validation or a similar repeated design.
#' predictionResults.txt contains a (data.frame)
#' @param inDir (char) path to predictionResults.txt files.
#' if inDir is a single char and not a vector of prediction files, the expected
#' directory structure is <inDir>/rng<X>/predictionResults.txt".
#' Otherwise inDir is assumed
#' to be a vector, each with absolute paths to predictionResults.txt
#' @param predClasses (char) vector of class names.
#' @return (list) each key corresponds to an input file in inDir.
#' Value is a list with:
#' 1) stats: "stats" component of perfCalc
#' 2) rocCurve: ROCR performance object for ROC curve
#' 3) prCurve: ROCR performance object for PR curve
#' 4) auroc: Area under ROC curve
#' 5) aupr: Area under PR curve
#' 6) accuracy: Accuracy
#'
#' Side effect of plotting in a 2x2 format:
#' 1) mean+/-SEM AUROC
#' 2) mean+/-SEM AUPR
#' 3) ROC curve for all runs plus average
#' 4) PR curve for all runs plus average
#' @examples
#' inDir <- sprintf("%s/extdata/KIRC_output",
#' 	path.package("netDx.examples"))
#' plotPerf(inDir, predClasses = c("SURVIVEYES", "SURVIVENO"))
#' @import ROCR
#' @import pracma
#' @export
plotPerf <- function(inDir, predClasses) {

  if (missing(inDir)) stop("inDir not provided");
  if (missing(predClasses)) stop("predClasses missing; please specify classes");

  # given output of performance("precall") compute AUC-PR
  prauc <- function(dat) {
  	x <- dat@x.values[[1]] # recall
  	y <- dat@y.values[[1]] # precision

  	# remove NAN
  	idx <- which(is.nan(y))
  	if (any(idx)) { x <- x[-idx]; y <- y[-idx]}

  	pracma::trapz(x,y)
  }

  if(length(inDir) == 1){
    cat("Single directory provided, retrieving prediction files\n")
    all_rng <- list.files(path = inDir, pattern = "rng.")
		fList <- sprintf("%s/%s/predictionResults.txt", inDir,all_rng)
	} else {
		cat("length(inDir)>1; assuming provided path to individual results\n")
		fList <- inDir
	}

  mega <- list()
  for (fName in fList) {
      out <- list()
        overall_acc <- numeric()
        curRoc	<- list()
        curPr	<- list()

        dat <- read.delim(fName,sep="\t",header=TRUE,as.is=TRUE)

        pred_col1 <- sprintf("%s_SCORE",predClasses[1])
        pred_col2 <- sprintf("%s_SCORE",predClasses[2])

			idx1 <- which(colnames(dat) == pred_col1)
			idx2 <- which(colnames(dat) == pred_col2)
        pred <- ROCR::prediction(dat[,idx1]-dat[,idx2],
					dat$STATUS==predClasses[1])

        c1 <- predClasses[1] #numc[1]
        tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
        tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
        fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
        fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)

		  # entire curves
        curRoc <- ROCR::performance(pred,"tpr","fpr")
        curPr <- ROCR::performance(pred,"prec","rec")
        tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn)
        out <- perfCalc(tmp)

        # statistic
        auroc	<- performance(pred, "auc")@y.values[[1]]
        aupr 	<- prauc(curPr)
        overall_acc <- c(overall_acc,
             sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)

### TODO put in F1.
      mega[[fName]] <- list(stats=out$stats,roc_curve=curRoc,pr_curve=curPr,
								auroc=auroc, aupr=aupr,accuracy=overall_acc)
    }
		# ----------------------------------
	  # Plot mean+/- SEM

		.plotAvg <- function(res, name) {
		 	mu <- mean(res,na.rm=TRUE)
  		sem <- sd(res,na.rm=TRUE)/sqrt(length(res))
 			plot(1, mu,type='n',bty='n',ylab=sprintf("%s (mean+/-SEM)",name),
       	 xaxt='n',ylim=c(0.4,1.0),las=1,cex.axis=1.4,xlim=c(0.8,1.2),
				cex.axis=1.4,xlab="")
    	abline(h=c(0.7,0.8),col='cadetblue3',lty=3,lwd=3)
    	points(1,mu,type='p',cex=1.4,pch=16)

    	# error bars
    	segments(x0=1, y0=mu-sem,y1=mu+sem,lwd=3)
    	segments(x0=1-0.01, x1=1+0.01,y0=mu-sem,y1=mu-sem)
    	segments(x0=1-0.01, x1=1+0.01,y0=mu+sem,y1=mu+sem)
    	abline(h=0.5,col='red',lty=1,lwd=2)
    	title(sprintf("%s: N=%i runs",name,length(res)))
		}

		# plot average +/-SEM
		
		par(mfrow=c(2,2))
		x <- unlist(lapply(mega, function(x) x$auroc))
		.plotAvg(x,"AUROC")
		x <- unlist(lapply(mega, function(x) x$aupr))
		.plotAvg(x,"AUPR")

		# plot individual curves
		rocCurves <- lapply(mega, function(x) x$roc_curve)
		plotPerf_multi(rocCurves,"ROC")
		prCurves <- lapply(mega, function(x) x$pr_curve)
		plotPerf_multi(prCurves,"PR",plotType="PR")

		return(mega)
}
