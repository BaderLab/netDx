#' Run a GeneMANIA query
#'
#' @param dbPath (char) path to directory with GeneMANIA generic database
#' @param queryFiles (list(char)) paths to query files
#' @param resDir (char) path to output directory
#' @param verbose (logical) print messages
#' @param JavaMemory (integer) Memory for GeneMANIA (in Gb) - a total of 
#' numCores*GMmemory will be used and distributed for all GM threads
#' @param numCores (integer) number of CPU cores for parallel processing
#' @return path to GeneMANIA query result file
#' of results file
#' @examples
#' dbPath <- sprintf("%s/extdata/dbPath", path.package("netDx"))
#' GM_query <- sprintf("%s/extdata/GM_query.txt",
#'		path.package("netDx"))
#' runQuery(dbPath, GM_query,"/tmp")
#' @export
runQuery <- function(dbPath, queryFiles, resDir, verbose=TRUE,
	JavaMemory=6L, numCores=1) {
	GM_jar	<- sprintf("%s/java/genemania-cytoscape-plugin-3.5.0.jar", path.package("netDx"))
	qBase	<- basename(queryFiles[[1]][1])
	logFile	<- sprintf("%s/%s.log", resDir, qBase)
	cmd1	<- sprintf("java -d64 -Xmx%iG -cp %s org.genemania.plugin.apps.QueryRunner",JavaMemory*numCores,GM_jar)
	queryStrings <- paste(queryFiles, collapse = ' ')
	cmd2	<- sprintf(" --data %s --in flat --out flat --threads %i --results %s %s --netdx-flag true 2>&1 > %s",
			dbPath, numCores, resDir, queryStrings, logFile)

	cmd		<- paste(c(cmd1,cmd2),collapse=" ")
	print(cmd)

	  # file is not actually created - is already split in PRANK and NRANK segments on GeneMANIA side
	resFile <- sprintf("%s/%s-results.report.txt", resDir,qBase)
	attempt <- 1
	# sometimes GM stochastically fails because of a failure-to-acquire-lock
	# in the /dataset/user directory. SP attributes to a race condition
	# when GM is executed in parallel. However, the problem did not seem
	# to occur in Feb 2016 on the VM but does occur in April 2016. 
	# This while-loop exists to ensure that all GM queries run.
	# removed in version 1.1
	t0	<- Sys.time()
	system(cmd,wait=TRUE,ignore.stdout=!verbose, ignore.stderr=!verbose)
	cat(sprintf("QueryRunner time taken: %1.1f s\n", Sys.time()-t0))

	Sys.sleep(1)
	return(resFile)
}
