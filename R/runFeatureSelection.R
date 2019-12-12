#' Run GeneMANIA cross-validation with a provided subset of networks
#'
#' @details Creates query files, runs GM for 10-fold cross validation.
#' @param trainID_pred (char) vector with universe of predictor class
#' patients (ie all that can possibly be included in the query file
#' @param outDir (char) directory to store query file and GM results
#' @param dbPath (char) path to GeneMANIA generic database with
#'\ttraining population
#' @param numTrainSamps (integer) number of training samples in total
#' leave blank to use 5 training samples in order to save memory
#' @param incNets (char) vector of networks to include in this analysis
#' (features/pathway names). Useful for subset-based feature selection
#' @param orgName (char) organism name for GeneMANIA generic database.
#' The default value will likely never need to be changed.
#' @param fileSfx (char) file suffix
#' @param verbose (logical) print messages
#' @param numCores (logical) num parallel threads for cross-validation
#' @param JavaMemory (integer) memory for GeneMANIA run, in Gb.
#' @param verbose_runQuery (logical) print messages for runQuery()
#' @param ... args for \code{makeQueries()}
#' @return No value. Side effect of generating feature scores.
#' @examples
#' data(MB.pheno)
#' dbPath <- sprintf('%s/extdata/dbPath',path.package('netDx'))
#' runFeatureSelection(MB.pheno$ID[which(MB.pheno$STATUS%in% 'WNT')],
#'\t'~/tmp',dbPath,103L)
#' @export
runFeatureSelection <- function(trainID_pred, outDir, dbPath, 
		numTrainSamps = NULL, incNets = "all", orgName = "predictor", 
		fileSfx = "CV", verbose = FALSE, numCores = 2L, 
    JavaMemory = 6L, verbose_runQuery = FALSE, ...) {
    
    if (!file.exists(outDir)) 
        dir.create(outDir)
    
    # get query names
    if (verbose) 
        message("\tWriting queries:\n")
    qSamps <- makeQueries(trainID_pred, verbose = verbose, ...)
    
    # write query files
    for (m in seq_len(length(qSamps))) {
        qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx, m)
        
        if (is.null(numTrainSamps)) {
            numTrainSamps = 5
            message("Memory saver option: using 5 training samples for CV")
        }
        
        writeQueryFile(qSamps[[m]], incNets, numTrainSamps, qFile, orgName)
    }
    qFiles <- list()
    for (m in seq_len(length(qSamps))) {
        qFile <- sprintf("%s/%s_%i.query", outDir, fileSfx, m)
        qFiles <- append(qFiles, qFile)
    }
    
    runQuery(dbPath, qFiles, outDir, JavaMemory = JavaMemory, 
				verbose = verbose_runQuery, 
        numCores = numCores)
    
}
