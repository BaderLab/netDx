#' Counts the relative correlation of (+,+) and (+,-)(-,-) interactions
#' 
#' @param inFile (character): path to interaction networks
#' @param plusID (character) vector of + nodes
#' @param minusID (character) vector of - nodes
#' @return (numeric) mean edge weight for (+,+) and other edges
getCorrType <- function(inFile, plusID, minusID) { 
	dat		<- read.delim(inFile,sep="\t",header=F,as.is=T)
	pp_idx	<- dat[,1] %in% plusID & dat[,2] %in% plusID
	pp_corr <- mean(dat[pp_idx,3])
	pm_corr <- mean(dat[setdiff(1:nrow(dat),pp_idx),3])

	return(c(pp_corr,pm_corr))
}
