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
#' @return (list) one key per patient class. Value is matrix of network
#' scores across all train/test splits. Each score is the output of
#' the inner fold of CV.
#' @examples
#' inDir <- sprintf("%s/extdata/KIRC_output",
#' 		path.package("netDx.examples"))
#' netScores <- getFeatureScores(inDir, predClasses = c("SURVIVEYES","SURVIVENO"))
#' @export
getFeatureScores <- function(inDir,predClasses) {
	if (missing(inDir)) stop("inDir not provided");
	if (missing(predClasses))
		stop("predClasses missing; please specify classes");

	out <- list()
	for (gp in predClasses) {
		cat(sprintf("%s\n",gp))

		if(class(inDir) == "character"){
			cat("\tSingle directory provided, retrieving CV score files\n")
			rngDirs <- dir(path=inDir, pattern="^rng")
			fList <-sprintf("%s/%s/%s/GM_results/%s_pathway_CV_score.txt",
					 inDir,rngDirs,gp,gp)
		} else {
			cat("\tList of filenames provided\n")
			fList <- inDir[[gp]]
		}

		cat(sprintf("Got %i iterations\n", length(fList)))
		netColl <- list()

		for (scoreFile in fList) {
			tmp	 <- read.delim(scoreFile,sep="\t",h=T,as.is=T)
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
			cat("* Computing consensus\n")
			cons <- getNetConsensus(netColl); x1 <- nrow(cons)
			na_sum <- rowSums(is.na(cons))
			full_cons <- cons
			cons <- cons[which(na_sum < 1),]

			out[[gp]] <- cons
	}
	return(out)
}
