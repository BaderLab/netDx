#' Feature subset selection driver function
#' 
#' @details Wrapper to run feature subselection methods from the FSelector 
#' package
#' @param methodName (char) subset selection strategy. Options are:
#' 1) greedy.forward : forward selection
#' 2) greedy.backward : backward elimination
#' 3) hill.climb: hill climbinb
#' @param incNets (char) vector of networks that are the basis of subset
#' selection
#' @param queryPool (char) See GM_run_featureSet, trainID param
#' @param GM_db (char) path to GeneMANIA generic database with 
#'	training population. This is the same as the \code{dbDir} value for
#' the output of compileFeatures.R
#' @param pheno_DF (data.frame) pheno matrix; ID and STATUS
#' @param predictClass (char) class label for predictor
#' @param outDir (char) directory to store query file and GM results
#' @param num2return (integer) number of training samples in total
#' @param numCores (integer) num cores for parallel processing. Used by
#' \code{runFeatureSelection()} during the cross-validation for each 
#' feature subset
#' @param scoreFile (char) path to file where error is printed over 
#' iterations. useful for debugging. If NULL doesn't print.
#' @export
#' @return (char) vector of selected subset of networks
#' @examples
#' GM_db <- sprintf("%s/extdata/GM_db",path.package("netDx"))
#' data(MB_pheno)
#' nets <- c("DKK2_cont","EMX2_cont","TNC_cont","WIF1_cont")
#' doSubsetSelection(methodName="greedy.backward",incNets=nets, 
#'	queryPool=MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
#'	GM_db=GM_db,pheno_DF=MB.pheno, outDir=".",num2return=103L)
doSubsetSelection <- function(methodName="greedy.forward",incNets,queryPool,
	GM_db,pheno_DF,predictClass,outDir,num2return,numCores=1L,
	scoreFile=NULL) {

	if (!is.null(scoreFile)) {
		cat("Feat score\n", file=scoreFile)
	}
	cat("Started at:\n"); print(Sys.time());

	# function run for each round for feature evaluation
	featEval <- function(netSubset) {
		t1 <- Sys.time();
		cat(sprintf("%i features\n",length(netSubset)))
		#curFile <- sprintf("%s/currently_eval_features.txt",outDir)
		#write.table(netSubset,file=curFile,sep="\t",
		#			row=F,col=F,quote=F)
	
		resDir <- sprintf("%s/subsetWork", outDir)
		runFeatureSelection(queryPool, resDir, GM_db, 
				num2return,
				verbose=FALSE,numCores=numCores,
				incNets=netSubset)

		prlist <- c()
		fList <- paste(resDir,dir(path=resDir,pattern="PRANK$"),sep="/")

		# TODO speed up by implementing precall compute
		# in here, instead of calling getPatientRankings.
		# can also parallelize
		for (f in fList)  {
			tmp <- getPatientRankings(f,pheno_DF,predictClass)
			tmp	<- tmp$precall
			df  <- data.frame(x=tmp@x.values[[1]],y=tmp@y.values[[1]])
			df	<- na.omit(df)
			prlist <- c(prlist,pracma::trapz(df$x,df$y))
		}
		err <- mean(prlist)

		if (!is.null(scoreFile))
			cat(sprintf("%1.3f\n", err),file=scoreFile,append=TRUE)
		t2 <- Sys.time()
		print(t2-t1)
		err
	 }

	fsub <- switch(methodName,
		greedy.forward={
			cat("* Running greedy forward selection\n")
			FSelector::forward.search(incNets,featEval)
		}, greedy.backward={
			cat("* Running greedy backward elimination\n")
		warning("not tested")
			FSelector::backward.search(incNets,featEval)
		}, hill.climb={
			cat("* Running hill climbing\n")
			FSelector::hill.climbing.search(incNets,featEval)
		}, stop("Invalid methodName. Type ?doSubsetSelection to see options"))
	cat("Finished at:\n"); print(Sys.time())
	fsub
}
