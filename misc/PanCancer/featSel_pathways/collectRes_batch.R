#' collect AUCROC and AUCPR for all conditions pathway related
require(ROCR)

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

# gets AUCROC and AUCPR from a table
.getAUC <- function(curFile) {
dat <- read.delim(curFile,sep="\t",h=T,as.is=T)
aucroc <- NA; aucpr <- NA
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
	aucroc <- performance(pred, "auc")@y.values[[1]]
	# precision-recall
	tmp <- performance(pred,"prec","rec")
	aucpr <- prauc(tmp)
}
return(c(aucroc,aucpr))
}

# --------------------------------------------------------------
# list of all conditions to collect data for and their i/o locations
setFile <- "KIRCpathway_locations.txt"
rootDir <- "/Users/shraddhapai/DropBox/netDx/BaderLab"
inRoot <- sprintf("%s/2017_TCGA_KIRC/output",rootDir)
outRoot <- sprintf("%s/2017_PanCancer_Survival",rootDir)

# --------------------------------------------------------------
dt <- format(Sys.Date(),"%y%m%d")
setInfo <- read.delim(setFile,sep="\t",h=T,as.is=T)
setInfo <- subset(setInfo, inc=="yes")

outList <- list()
megaList <- list()

for (ctr in 1:nrow(setInfo)) {
	cur <- setInfo$name[ctr]
	curd <- sprintf("%s/%s",inRoot,setInfo$dataDir[ctr])
	maxk <- setInfo$maxK[ctr]
	outDir <- sprintf("%s/%s", outRoot, setInfo$outdir[ctr])
	if (!file.exists(outDir)) dir.create(outDir)
	outFile <- sprintf("%s/%s/KIRC_results_%s.Rdata",
		outRoot,setInfo$outdir[ctr],cur)

	kset <- 1:maxk
	cat(sprintf("%s:Num runs=%i\n", cur,maxk))

	val		<- matrix(NA,nrow=maxk,ncol=1)
	val_pr <- matrix(NA,nrow=maxk,ncol=1)
	colnames(val) <- ""
	if (cur == "pathOnlyRnd") { # random nets have different structure
		dirs <- dir(curd,pattern="random")
		big_ctr <- 1
		for (d in dirs) { # one round of sampling
			fSet <- dir(sprintf("%s/%s/predictions",curd,d),pattern="prediction")
			cat(sprintf("Got %i predictions\n", length(fSet)))
			tmp <- matrix(NA,nrow=length(fSet),ncol=2)
			ctr <- 1
			for (fName in fSet) {	
					curF <- sprintf("%s/%s/predictions/%s",curd,d,fName)
					tmp[ctr,] <- .getAUC(curF)
					ctr <- ctr+1
			}
			val[big_ctr,1] <- mean(tmp[,1])
			val_pr[big_ctr,1] <- mean(tmp[,2])
			big_ctr <- big_ctr + 1
		}
} else {
	for (k in kset) {
			print(k)
			if (cur == "pathOnlyCons") {
				curFile	<-sprintf("%s/predictions/predictionResults_%i.txt",curd,k)
			} else { 
				curFile<-sprintf("%s/rng%i/predictionResults.txt",curd,k)
			}
			tmp <- .getAUC(curFile);
			val[k,1] <- tmp[1]
			val_pr[k,1] <- tmp[2]
	}
}

	cat("Saving to file.\n")
	save(val,val_pr,file=outFile)
}
