#' Convert profiles to interaction networks before integration
#'
#' @details In preparation for network integration. When using GeneMANIA's
#' built-in functionality to create PSN using ProfileToNetworkDriver, this
#' step needs to run to process profiles to networks. These are currently used
#' for Pearson correlation-based networks and those using mutual information.
#' @param netDir (char) directory with .profile files
#' @param outDir (char) path to directory where interaction networks are to be printed
#' @param simMetric (char) similarity measure to use in converting 
#' profiles to interaction networks. 
#' @param numCores (integer) number of cores for parallel processing 
#' @param P2N_threshType (char) Most users shouldn't have to change this.
#' ProfileToNetworkDriver's threshold option. One of 'off|auto'. 
#' unit testing
#' @param P2N_maxMissing (integer 5-100)
#' @param JavaMemory (integer) Memory for GeneMANIA (in Gb)
#' @param GM_jar (char) path to GeneMANIA jar file
#' @param netSfx (char) pattern for finding network files in \code{netDir}.
#' @param debugMode (logical) if TRUE runs profile generation in serial 
#' rather than parallel, allowing debugging
#' @return No value. Side effect of creating interaction networks in outDir.
#' @export
convertProfileToNetworks <- function(netDir,outDir=tempdir(),
	simMetric="pearson",numCores=1L,
	JavaMemory=4L,GM_jar=NULL,P2N_threshType="off",P2N_maxMissing=100,
	netSfx="txt$",debugMode=FALSE) {

if (is.null(GM_jar)) GM_jar <- getGMjar_path()

cl <- makeCluster(numCores, 
	outfile = paste(netDir, "P2N_log.txt",sep=getFileSep()))
registerDoParallel(cl)

if (simMetric == "pearson") {
	corType <- "PEARSON"
} else if (simMetric == "MI") {
	corType <- "MUTUAL_INFORMATION"
}
        
args <- c(sprintf("-Xmx%iG", JavaMemory), "-cp", GM_jar)
args <- c(args, 
	paste("org.genemania.engine.core.",
	"evaluation.ProfileToNetworkDriver",sep=""))
args <- c(args, c("-proftype", "continuous", "-cor", corType))
args <- c(args, c("-threshold", P2N_threshType, 
	"-maxmissing", 
	sprintf("%1.1f", P2N_maxMissing)))
profDir <- netDir
tmpsfx <- sub("\\$", "", netSfx)
        
curProf <- ""

`%myinfix%` <- ifelse(debugMode, `%do%`, `%dopar%`)
foreach(curProf = dir(path = profDir, pattern = "profile$")) %myinfix% {
	if (debugMode) print(curProf)
	args2 <- c("-in", paste(profDir, curProf,sep=getFileSep()))
	args2 <- c(args2, "-out", paste(outDir, 
		sub(".profile", ".txt", curProf),sep=getFileSep()))
	args2 <- c(args2, "-syn", 
		paste(netDir,"..","1.synonyms",sep=getFileSep()),
		"-keepAllTies", "-limitTies")

	if (debugMode) stdout <- "" else stdout <- NULL
	system2("java", args = c(args, args2), wait = TRUE, 
		stdout = stdout)
}
stopCluster(cl)

}
