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
runQuery <- function(dbPath, queryFiles, resDir, verbose = TRUE,
    JavaMemory = 6L, numCores = 1L, debugMode = FALSE) {

  GM_jar <- getGMjar_path()
  qBase <- basename(queryFiles[[1]][1])
  logFile <- paste(resDir, sprintf("%s.log", qBase))
  queryStrings <- paste(queryFiles, collapse = " ")

  args <- c()
  java_ver <- suppressMessages(suppressWarnings(system2("java",
    args = "--version", stdout = TRUE, stderr = NULL)))
  if (any(grep(" 11", java_ver)) ||
      any(grep(" 12", java_ver)) ||
      any(grep(" 13", java_ver)) ||
      any(grep(" 14", java_ver)) ||
      any(grep(" 16", java_ver))) ||
      any(grep(" 17", java_ver)) ||
      any(grep(" 18", java_ver)) ||
	    any(grep(" 19", java_ver)) ||
	    any(grep(" 20", java_ver)) {
    if (verbose) message("Java 11 or later detected")
    args <- c("--illegal-access=permit") # needed for Java 9-16. Deprecated in Java 17)
  } else {
    if (verbose) message("Java 8 detected")
    args <- c(args, "-d64")
  }

  args <- c(args, sprintf("-Xmx%iG", JavaMemory * numCores), "-cp", GM_jar)
  args <- c(args, "org.genemania.plugin.apps.QueryRunner")
  args <- c(args, "--data", dbPath, "--in", "flat", "--out", "flat")
  args <- c(args, "--threads", numCores, "--results", resDir,
      unlist(queryFiles))
  args <- c(args, "--netdx-flag", "true") #,'2>1','/dev/null')

  # file is not actually created - is already split in PRANK and 
  # NRANK segments on
  # GeneMANIA side
  resFile <- paste(resDir, sprintf("%s-results.report.txt", qBase),
    sep = getFileSep())
  t0 <- Sys.time()
  if (debugMode) {
    message(sprintf("java %s", paste(args, collapse = " ")))
    system2("java", args, wait = TRUE)
  } else {
    system2("java", args, wait = TRUE, stdout = NULL, stderr = NULL)
  }
  if (verbose)
    message(sprintf("QueryRunner time taken: %1.1f s", Sys.time() - t0))
  Sys.sleep(3)
  return(resFile)
}
