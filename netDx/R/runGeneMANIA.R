#' Run a GeneMANIA query
#'
#' @param GM_db (char) path to directory with GeneMANIA generic database
#' @param queryFile (char) path to query file
#' @param resDir (char) path to output directory
#' @param parseReport (logical) if TRUE, parses out PRANK and NRANK portions
#' @param verbose (logical) print messages
#' @param GMmemory (integer) Memory for GeneMANIA (in Gb)
#' @param MAX_ATTEMPTS (integer) max num attempts to run GeneMANIA before
#' giving up.
#' @return path to GeneMANIA query result file
#' of results file
#' @export
runGeneMANIA <- function(GM_db, queryFile, resDir, parseReport=TRUE,
	verbose=TRUE,GMmemory=6L,MAX_ATTEMPTS=3L) {
	GM_jar	<- sprintf("%s/java/GeneMANIA-3.2B7.jar",
						 path.package("netDx"))
	qBase	<- basename(queryFile)
	logFile	<- sprintf("%s/%s.log", resDir, qBase)
	cmd1	<- sprintf("java -d64 -Xmx%iG -cp %s org.genemania.plugin.apps.QueryRunner",GMmemory,GM_jar)
	cmd2	<- sprintf(" --data %s --in flat --out flat --threads %i --results %s %s 2>&1 > %s",
			GM_db, 1, resDir, queryFile,logFile)

	cmd		<- paste(c(cmd1,cmd2),collapse=" ")
	print(cmd)
	
	resFile <- sprintf("%s/%s-results.report.txt", resDir,qBase)
	attempt <- 1
	# sometimes GM stochastically fails because of a failure-to-acquire-lock
	# in the /dataset/user directory. SP attributes to a race condition
	# when GM is executed in parallel. However, the problem did not seem
	# to occur in Feb 2016 on the VM but does occur in April 2016. 
	# however, this while-loop exists to ensure that all GM queries run.
	while ((!file.exists(resFile)) & (attempt <= MAX_ATTEMPTS)) {
			cat(sprintf("* Attempt %i : %s\n", attempt,
						basename(queryFile)))
		system(cmd,wait=TRUE,ignore.stdout=!verbose, ignore.stderr=!verbose)
		attempt <- attempt + 1
	}
	
	Sys.sleep(3)
	if (parseReport) GM_parseReport(resFile)
	
	return(resFile)
}
