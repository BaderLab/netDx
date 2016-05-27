#' Counts the number of (+,+) and (+,-) interactions in a single network
#' 
#' @param inFile (char) path to interaction networks
#' @param plusID (char) vector of + nodes
#' @param minusID (char) vector of - nodes
countIntType <- function(inFile, plusID, minusID) { 
	dat <- read.delim(inFile,sep="\t",header=F,as.is=T)
	pp	<- sum(dat[,1] %in% plusID & dat[,2] %in% plusID)

	return(c(pp,nrow(dat)-pp))
}
