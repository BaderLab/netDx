#' Counts the number of (+,+) and (+,-) interactions in a single network
#' 
#' @param inFile (char) path to interaction networks
#' @param plusID (char) vector of + nodes
#' @param minusID (char) vector of - nodes
#' @return (numeric of length 2) Number of (+,+) interactions, and non-(+,+) interactions
#' (i.e. (+,-) and (-,-) interactions)
#' @export
#' @examples
#' data(npheno)
#' netDir <- sprintf("%s/extdata/example_nets",
#'		path.package("netDx"))
#' countIntType(sprintf("%s/BOTH_EQUAL.txt", netDir),
#' 		npheno[1:100,1],npheno[101:200,1])
countIntType <- function(inFile, plusID, minusID) { 
	dat <- read.delim(inFile,sep="\t",header=F,as.is=T)
	pp	<- sum(dat[,1] %in% plusID & dat[,2] %in% plusID)

	return(c(pp,nrow(dat)-pp))
}
