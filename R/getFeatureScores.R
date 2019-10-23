#' Compile network scores into a matrix
#'
#' @details Given network scores over a set of trian/test splits, compiles these
#' into a matrix for downstream analysis. See the section on "Output Files"
#' @param inDir (char/list) directory containing directories with all split info
#' or list of all CV score files.
#' if inDir is a single directory then the expected format for CV score files is
#' <inDir>/rngX/predClassX/GM_results/predClassX_pathway_CV_score.txt"
#' if inDir is a list, it should have one key per class. The value should be the
#' corresponding set of filenames for pathway_CV_score.txt
#' @param predClasses (char) possible STATUS for patients
#' @param getFullCons (logical) if TRUE, does not remove rows with NA.
#' Recommended only when the number of input features is extensively 
#' pruned by first-pass feature selection.
#' @return (list) one key per patient class. Value is matrix of network
#' scores across all train/test splits. Each score is the output of
#' the inner fold of CV.
#' @examples
#' inDir <- sprintf("%s/extdata/example_output",path.package("netDx"))
#' netScores <- getFeatureScores(inDir, predClasses = c("SURVIVEYES","SURVIVENO"))
#' @export
getFeatureScores <- function(inDir,predClasses,getFullCons=FALSE) {
	if (missing(inDir)) stop("inDir not provided");
	if (missing(predClasses))
		stop("predClasses missing; please specify classes");

	out <- list()
	for (gp in predClasses) {
		message(sprintf("%s\n",gp))

		if(is(inDir,"character")) {
			message("\tSingle directory provided, retrieving CV score files\n")
			rngDirs <- dir(path=inDir, pattern="^rng")
			fList <-sprintf("%s/%s/%s/GM_results/%s_pathway_CV_score.txt",
					 inDir,rngDirs,gp,gp)
		} else {
			message("\tList of filenames provided\n")
			fList <- inDir[[gp]]
		}

		message(sprintf("Got %i iterations\n", length(fList)))
		netColl <- list()

		for (scoreFile in fList) {
			tmp	 <- read.delim(scoreFile,sep="\t",header=TRUE,as.is=TRUE)
			colnames(tmp)[1] <- "PATHWAY_NAME"
				netColl[[scoreFile]] <- tmp
		}
			spos <- gregexpr("\\/",fList)
			# get the name of the iteration (rngX) assuming directory structure
			# rngX/<class>/GM_results>/pathway_CV_score.txt
			fNames <- lapply(1:length(spos), function(x) {
				  n <- length(spos[[x]])
					y <- substr(fList[x], spos[[x]][n-3]+1,spos[[x]][n-2]-1)
					y
			})
			fNames <- unlist(fNames)
			names(netColl) <- fNames

			# filter for nets meeting cutoff criteria
			message("* Computing consensus\n")

			cons <- getNetConsensus(netColl); x1 <- nrow(cons)
			na_sum <- rowSums(is.na(cons))
			full_cons <- cons
			cons <- cons[which(na_sum < 1),]

			if (getFullCons) out[[gp]] <- full_cons else out[[gp]] <- cons
	}
	return(out)
}
