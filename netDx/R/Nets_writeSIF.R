#' write patient networks in Cytoscape's .sif format
#'
#' @details Converts a set of binary interaction networks into Cytoscape's
#' sif format.
#' (http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats)
#' This utility permits visualization of feature selected networks.
#'
#' @param netPath (char): vector of path to network files; file suffix
#' should be "_cont.txt" 
#' networks should be in format:	A	B	1
#' where A and B are nodes, and 1 indicates an edge between them
#' @param outFile (char) path to .sif file 
#' @param netSfx (char) suffix for network file name
#' @return No value. Side effect of writing all networks to \code{outFile}
#' @examples
#' netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
#' netFiles <- sprintf("%s/%s", netDir, dir(netDir,pattern="txt$"))
#' writeNetsSIF(netFiles,"merged.sif",netSfx=".txt")
#' @export
writeNetsSIF <- function(netPath,outFile,netSfx="_cont.txt"){
	
	system(sprintf("cat /dev/null > %s",outFile))
	for (n in netPath) {
		netName <- sub(netSfx,"",basename(n))
		cat(sprintf("%s\n", netName))

		dat <- read.delim(n,sep="\t",h=F,as.is=T)
		dat2 <- cbind(dat[,1],netName,dat[,2])

		write.table(dat2,file=outFile,append=TRUE,sep="\t",col=F,row=F,
					quote=F)
	}
	
}
