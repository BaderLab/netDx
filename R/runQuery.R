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
	JavaMemory=6L, numCores=1L) {
	
	GM_jar	<- getGMjar_path()
	qBase	<- basename(queryFiles[[1]][1])
	logFile	<- sprintf("%s/%s.log", resDir, qBase)
	queryStrings <- paste(queryFiles, collapse = ' ')
	#cmd1	<- sprintf("java -d64 -Xmx%iG -cp %s org.genemania.plugin.apps.QueryRunner",JavaMemory*numCores,GM_jar)
	#cmd2	<- sprintf(" --data %s --in flat --out flat --threads %i --results %s %s --netdx-flag true",
	#		dbPath, numCores, resDir, queryStrings)
	#cmd		<- paste(c(cmd1,cmd2),collapse=" ")
	#if (!verbose) cmd <- sprintf("%s 2>1 /dev/null",cmd)
	#if (verbose) print(cmd)
	#blah <- suppressWarnings(system2(cmd,wait=TRUE,stdout=TRUE))

	args <- c('-d64',sprintf('-Xmx%iG',JavaMemory*numCores),'-cp',GM_jar)
	args <- c(args,'org.genemania.plugin.apps.QueryRunner')
	args <- c(args,'--data', dbPath, '--in','flat','--out','flat')
	args <- c(args,'--threads',numCores,'--results',resDir,unlist(queryFiles))
	args <- c(args,'--netdx-flag','true') #,'2>1','/dev/null')

	# file is not actually created - is already split in PRANK and NRANK 
	# segments on GeneMANIA side
	resFile <- sprintf("%s/%s-results.report.txt", resDir,qBase)
	t0	<- Sys.time()
	system2('java',args,wait=TRUE,stdout=NULL,stderr=NULL)
	if (verbose) message(sprintf("QueryRunner time taken: %1.1f s", 
		Sys.time()-t0))
	Sys.sleep(3)
	return(resFile)
}
