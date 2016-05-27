#' Run GeneMANIA cross-validation with a provided subset of networks
#'
#' @details Creates query files, runs GM for 10-fold cross validation. 
#' @param trainID_pred (char) vector with universe of predictor class 
#' patients (ie all that can possibly be included in the query file
#' @param outDir (char) directory to store query file and GM results
#' @param GM_db (char) path to GeneMANIA generic database with 
#'	training population
#' @param numTrainSamps (integer) number of training samples in total
#' @param incNets (char) vector of networks to include in this analysis 
#' (features/pathway names). Useful for subset-based feature selection
#' @param orgName (char) organism name for GeneMANIA generic database.
#' The default value will likely never need to be changed.
#' @param fileSfx (char) file suffix
#' @param verbose (logical) print messages
#' @param numCores (logical) num parallel threads for cross-validation
#' @param GMmemory (integer) memory for GeneMANIA run, in Gb.
#' @param ... args for \code{makeCVqueries()}
#' @export
GM_runCV_featureSet <- function(trainID_pred,outDir,GM_db,numTrainSamps, 
	incNets="all",orgName="predictor",fileSfx="CV",verbose=FALSE,
	numCores=2L,GMmemory=6L,...) {	
	
	#TODO if results already exist, what do we do? Delete with a warning?
	if (!file.exists(outDir)) dir.create(outDir)

	# get query names 
	if (verbose) cat("\tWriting GM queries: ")
	qSamps <- makeCVqueries(trainID_pred,verbose=verbose,...)

	# write query files
	for (m in 1:length(qSamps)) {
		if (verbose) cat(sprintf("%i ",m))
		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx,m)
		GM_writeQueryFile(qSamps[[m]], incNets, numTrainSamps, 
						  qFile,orgName)
	}

	cl	<- makeCluster(numCores)
	registerDoParallel(cl)

	# run GeneMANIA 10-fold
	x <- foreach(m=1:length(qSamps)) %dopar% {
		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx,m)

		runGeneMANIA(GM_db, qFile, outDir,GMmemory=GMmemory,
					 verbose=verbose)
		
		# needed so R will pause while GM finishes running. 
		# otherwise the script proceeds asynchronously (without waiting
		# for GM to complete) and will try to work with result files that
		# haven't been written yet
		#Sys.sleep(15) 	
	
		# keep only PRANK and NRANK and remove main results file.
		resFile <- sprintf("%s/%s_%i.query-results.report.txt",
			outDir,fileSfx,m)

		system(sprintf("rm %s",resFile))
	}

	stopCluster(cl)
}
