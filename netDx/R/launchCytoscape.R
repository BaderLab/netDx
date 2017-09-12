#' launches Cytoscape
#'
#' @details Known issues. The output associated with Cytoscape launch causes
#' the terminal to be useless for further interaction because although it
#' executes typed commands, it doesn't print the characters to the screen.
#' This script should only be used in completely non-interactive sessions.
#' @param cytoPath (char) Absolute path to installed Cytoscape software
#' @param sleepTime (integer) number of seconds to sleep while waiting for
#' Cytoscape to launch. If cyrest tries to run commands too soon, you may
#' want to increase this duration.
#' @return No value. Side effect of launching Cytoscape.
#' @example
#' launchCytoscape()
#' @export
launchCytoscape <- function(cytoPath="/Applications",sleepTime=30) {

cat("Launching Cytoscape\n")
x <- dir(cytoPath,pattern="Cytoscape_")

if (length(x)<1) {
	cat("-------------------------------------------------------\n")
	cat("Cytoscape cannot be found in the install path provided.\n")
	cat(sprintf("Install path: %s\n", cytoPath))
	cat("If Cytoscape is installed on another path in your system,\n")
	cat("rerun this command by changing the cytoPath variable.\n")
	cat("\n")
	cat("If Cytoscape is not installed, visit http://cytoscape.org\n")
	cat("to install before continuing.\n")
	cat("-------------------------------------------------------\n")
	stop()
}

if (length(x)>1) {
	x <- x[length(x)]
	cat(sprintf("\tMultiple versions detected; picking %s",x))
}
x <- sprintf("%s/%s/cytoscape.sh",cytoPath, x)
system(sprintf("sh %s", x),wait=FALSE)
cat(sprintf("\n\nPausing %i seconds for Cytoscape to start up (one-time cost)...\n",
	sleepTime))
Sys.sleep(sleepTime) # wait for Cytoscape to start

}
