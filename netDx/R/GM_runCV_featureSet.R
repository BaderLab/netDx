#' Run GeneMANIA cross-validation with a provided subset of networks
#'
#' @details Creates query files, runs GM for 10-fold cross validation.
#' @param trainID_pred (char) vector with universe of predictor class
#' patients (ie all that can possibly be included in the query file
#' @param outDir (char) directory to store query file and GM results
#' @param GM_db (char) path to GeneMANIA generic database with
#'	training population
#' @param numTrainSamps (integer) number of training samples in total
#' leave blank to use 5 training samples in order to save memory
#' @param incNets (char) vector of networks to include in this analysis
#' (features/pathway names). Useful for subset-based feature selection
#' @param orgName (char) organism name for GeneMANIA generic database.
#' The default value will likely never need to be changed.
#' @param fileSfx (char) file suffix
#' @param verbose (logical) print messages
#' @param numCores (logical) num parallel threads for cross-validation
#' @param useNewGM (logical) use new GeneMania Version for processing
#' @param GMmemory (integer) memory for GeneMANIA run, in GB.
#' @param seed_CVqueries (integer) RNG seed for inner cross validation loop.
#' Makes deterministic samples held-out for each GeneMANIA query (see
#' makeCVqueries())
#' @param ... args for \code{makeCVqueries()}
#' @examples
#' data(MB_pheno)
#' GM_db <- sprintf("%s/extdata/GM_db",path.package("netDx"))
#' GM_runCV_featureSet(MB.pheno$ID[which(MB.pheno$STATUS%in% "WNT")],
#'	"~/tmp",GM_db,103L)
#' @export
GM_runCV_featureSet <- function(trainID_pred,outDir,GM_db,numTrainSamps = NULL,
	incNets="all",orgName="predictor",fileSfx="CV",verbose=FALSE,
	numCores=2L,
	useNewGM=FALSE,
	GMmemory=6L,seed_CVqueries=42L,...) {

	#TODO if results already exist, what do we do? Delete with a warning?
	if (!file.exists(outDir)) dir.create(outDir)

	# get query names
	if (verbose) cat("\tWriting GM queries: ")
	qSamps <- makeCVqueries(trainID_pred,verbose=verbose,
		setSeed=seed_CVqueries,...)

	# write query files
	for (m in 1:length(qSamps)) {
		if (verbose) cat(sprintf("%i ",m))
		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx,m)

		if(is.null(numTrainSamps)){
			numTrainSamps = 5
			cat("Memory saver option: using 5 training samples for CV")
		}

		GM_writeQueryFile(qSamps[[m]], incNets, numTrainSamps,
						  qFile,orgName)
	}
  if (useNewGM ){
    qFiles <- list()
    for (m in 1:length(qSamps)) {
      qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx, m)
      qFiles <- append(qFiles, qFile)
    }
    
    runGeneMANIA3(GM_db, qFiles, outDir,GMmemory=GMmemory,
                  verbose=verbose,parseReport=FALSE, numCores=numCores)
    
  } else{ 
	
  	cl	<- makeCluster(numCores)
  	registerDoParallel(cl)
  
  	# run GeneMANIA n-fold
  	x <- foreach(m=1:length(qSamps),
  	             .packages = c("netDx")) %dopar% {
  		qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx, m)
  
  		# its also possible to pass multiple / all query files at once to GeneMania 
  		# and increase the number of used threads in GeneMania --> useGMThreads
		  runGeneMANIA(GM_db, qFile, outDir,GMmemory=GMmemory,
  		               verbose=verbose)
  		  
  		  # keep only PRANK and NRANK and remove main results file.
  		  resFile <- sprintf("%s/%s_%i.query-results.report.txt",
  		                     outDir,fileSfx,m)
  		  
  		  system(sprintf("rm %s",resFile))
  		}
  		
  		# needed so R will pause while GM finishes running.
  		# otherwise the script proceeds asynchronously (without waiting
  		# for GM to complete) and will try to work with result files that
  		# haven't been written yet
  		#Sys.sleep(15)
  
  	stopCluster(cl)
  }
}
