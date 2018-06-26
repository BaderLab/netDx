#' compare the similarity ranks of test patients to each of the target
#' classes
#' Expects a directory structure like so:
#' <inDir>
#'  -- rng1/predictionResults.txt
#'  -- rng2/predictionResults.txt
#' ...
#' 
#' @param inDir (char) path to input directory
#' @param classes (char) patient classes
compareTestSimRanks <- function(inDir,classes=c("SURVIVEYES","SURVIVENO")) {
  dirSet <- dir(path=inDir,pattern="^rng")
	cat(sprintf("Got %i dirs\n",length(dirSet)))
	scoreList <- list()
	k <- 1
	for (d in dirSet) {
			#cat(sprintf("\t%s\n",d))
			dat <- read.delim(sprintf("%s/%s/predictionResults.txt",inDir,d),
				sep="\t",h=T,as.is=T)
			x1 <- dat[,paste(classes[1],"_SCORE",sep="")]
			x2 <- dat[paste(classes[2],"_SCORE",sep="")]
			df <- x1-x2
			scoreList[[k]] <- data.frame(ID=dat$ID, scoreDiff=df)
			colnames(scoreList[[k]])[2] <- "scoreDiff"
			k <- k+1
	}

	cat("Compiling single matrix\n")
	uq_pat <- unique(unlist(lapply(scoreList, function(x) x$ID)))
	mat <- matrix(NA,nrow=length(uq_pat),ncol=length(scoreList))
	rownames(mat) <- uq_pat
	
	for (k in 1:length(scoreList)) {
		cur <- scoreList[[k]]
		cur[,1] <- as.character(cur[,1])
		midx <- match(cur$ID,rownames(mat))
		if (all.equal(rownames(mat)[midx],cur$ID)!=TRUE) {
			cat("ids don't match\n")
			browser()
		}
		mat[midx,k] <- cur$scoreDiff
	}
	
	return(mat)
}
