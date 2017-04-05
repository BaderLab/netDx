# assess performance of BRCA classifier for different combinations of
# input data types
rm(list=ls())
require(netDx)

rootDir <- "/Users/shraddhapai/Documents/Research/BaderLab"
dataName <- "OV"
inDirList <- list(
	KIRC=sprintf("%s/2017_TCGA_KIRC/output/featSel_170222", rootDir),
	OV=sprintf("%s/2017_TCGA_OV/output/OV_170227",rootDir),
	LUSC=sprintf("%s/2017_TCGA_LUSC/output/featSel_incMutRPPA_round2170223",
							 rootDir),
	GBM=sprintf("%s/2017_TCGA_GBM/output/featSel_incMut_round2_170223",rootDir)
	)
inDir <- inDirList[[dataName]]

cat("--------------------\n")
cat(sprintf("%s\n",dataName))
cat("--------------------\n")

# combSet <- c("clinical","clinicalArna","clinicalAmir","clinicalArppa","all")
combSet <- c("all")
cols <- c(brewer.pal(n=4,name="Dark2"),"red")

# --------------------------------------------------------------
# helper routines

require(RColorBrewer)
.plotROC <- function(dat,pType) {
	plot(NA,NA,
		 ylab="",xlab="",xlim=c(0,1),ylim=c(0,1),bty='n')
	for (k in 1:length(dat)) {
		if (k==4) lwd=1.4 else lwd=1
		x <- dat[[k]]@x.values[[1]]
		y <- dat[[k]]@y.values[[1]]
		points(x,y,type='l',col=cols[k],lwd=lwd)
	}
	text(0.5,0,sprintf("N=%i",length(dat[[1]]@x.values[[1]])))

	box(col='grey80')
	if (pType == "roc") abline(0,1,col='grey50')
	else abline(h=0.5,col='grey50')
}

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
# begin work
cat(sprintf("Got %i combs\n", length(combSet)))

mega <- list()
for (rngSeed in setdiff(1:100,97)) {
	out <- list()
	overall_acc <- numeric()
	curRoc	<- list()
	curPr	<- list()
	for (cur in combSet) {
		# collect info for ROC curve
		#x <- outRes$same$roc; 
		#x2 <- outRes$orse$roc;
		#curRoc[[cur]] <- list(same=outRes$same$roc, 
		#			  worse=outRes$worse$roc)
		#curPr[[cur]] <- list(same=outRes$same$precall,
	#				 worse=outRes$worse$precall)

		inf <- sprintf("%s/rng%i/predictionResults.txt",
					   inDir,rngSeed)
		dat <- read.delim(inf,sep="\t",h=T,as.is=T)
		pred <- prediction(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE,
						  dat$STATUS=="SURVIVEYES")

		c1 <- "SURVIVEYES" #numc[1]
		tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
		tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
		fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
		fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)

		# dummy score column needed for perfCalc()
		tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn) 
		curRoc[[cur]] <- performance(pred,"tpr","fpr")
		curPr[[cur]] <- performance(pred,"prec","rec")

		out[[cur]] <- perfCalc(tmp)
		# add AUC-PR and AUCROC
		out[[cur]]$auc <- performance(pred, "auc")@y.values[[1]]
		out[[cur]]$prauc <- prauc(curPr[[cur]])
		#out[[cur]]$aucroc #<- list(same=outRes$same$auc, 
				#	  worse=outRes$worse$auc,
				#	  comb=mean(c(outRes$same$auc, 
				#	outRes$worse$auc)))
		#out[[cur]]$aucpr <- list(same=prauc(outRes$same$precall),
		#			 worse=prauc(outRes$worse$precall))
		#out[[cur]]$aucpr$comb <- mean(unlist(out[[cur]]$aucpr))
		overall_acc <- c(overall_acc, 
			 sum(dat$STATUS==dat$PRED_CLASS)/nrow(dat)*100)
	}
	names(overall_acc) <- combSet
	mega[[as.character(rngSeed)]] <- list(out=out,overall_acc=overall_acc,
				roc=curRoc,pr=curPr)
}

# get stat interest from data structure with all results
func <- function(str) {
	tmp <- lapply(mega, function(x) {
		y <- x$out
		z <- unlist(lapply(y,function(w) { w$stats[[str]]}))
		z
	})
	tmp <- do.call("rbind",tmp)

	# normalize to first column
	#for (k in ncol(tmp):1) tmp[,k] <- tmp[,k]-tmp[,1]
	tmp
}

# extract performance measures, compile to plotting format
aucroc	<- lapply(mega,function(x) {
	z <- lapply(x$out, function(y) y$auc)
	unlist(z)
	})
aucroc <- do.call("rbind",aucroc)

best_rng <- as.integer(order(aucroc)[length(aucroc)])
best_aucroc <- aucroc[best_rng]

cat(sprintf("best rng %i, AUCROC: %f\n", best_rng, best_aucroc))

blah <- unlist(lapply(mega, function(x) x$overall_acc))
cat(sprintf("Best Accuracy = %1.2f , round %i\n", max(blah),which.max(blah)))



