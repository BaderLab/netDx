#' Run a query
#'
#' @param dbPath (char) path to directory with GeneMANIA generic database
#' @param queryFiles (list(char)) paths to query files
#' @param resDir (char) path to output directory
#' @param verbose (logical) print messages
#' @param JavaMemory (integer) Memory for GeneMANIA (in Gb) - a total of 
#' numCores*GMmemory will be used and distributed for all GM threads
#' @param numCores (integer) number of CPU cores for parallel processing
#' @param debugMode (logical) when TRUE runs jobs in serial instead of parallel and 
#' prints verbose messages. Also prints system Java calls.
#' @return (char) path to GeneMANIA query result files with patient similarity
#' rankings (*PRANK) and feature weights (*NRANK)
#' of results file
#' @examples
#' dbPath <- system.file("extdata","dbPath",package="netDx")
#' queryFile <- system.file("extdata","GM_query.txt",package="netDx")
#' runQuery(dbPath, queryFile,tempdir())
#' @export
#' @importFrom rJava .jinit .jcheck .jaddClassPath .jcall .jnew .jclassPath
runQuery <- function(dbPath, queryFiles, resDir, verbose = TRUE, 
		JavaMemory = 6L, numCores = 1L,debugMode=FALSE) {
    
    GM_jar <- getGMjar_path()
    qBase <- basename(queryFiles[[1]][1])
    logFile <- paste(resDir,sprintf("%s.log",qBase))
    queryStrings <- paste(queryFiles, collapse = " ")

  java_ver <- .jcall("java/lang/System", "S",
    "getProperty", "java.runtime.version")
 dpos <- unlist(gregexpr("\\.",java_ver)[[1]])
 java_ver <- substr(java_ver, 1, dpos[2]-1)
    args <- c("--data", dbPath, "--in", 
	"flat", "--out", "flat")
    args <- c(args, "--threads", numCores, 
	"--results", resDir, 
	unlist(queryFiles))
    args <- c(args, "--netdx-flag", "true")  #,'2>1','/dev/null')
    
    # file is not actually created - is already split in PRANK and 
	# NRANK segments on
    # GeneMANIA side
    resFile <- paste(resDir,
	sprintf("%s-results.report.txt",qBase),
	sep=getFileSep())
    t0 <- Sys.time()
	x <- getGMjar_path()
   tryCatch({
        .jcheck(silent=FALSE)
    },error=function(ex){
        .jinit()
        .jaddClassPath(x)
    })
    if (!x %in% .jclassPath()) .jaddClassPath(x)
  jObj <- .jnew("org.genemania.plugin.apps.QueryRunner")
  .jcall(jObj,"V",method="main",args)

    if (verbose) 
        message(sprintf("QueryRunner time taken: %1.1f s", 
	Sys.time() - t0))
    Sys.sleep(3)
    return(resFile)
}
