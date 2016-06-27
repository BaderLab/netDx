#' resolve patient number in GeneMANIA interaction net to names
#'
#' @param netFile (char) path to GeneMANIA interaction net. SOURCE TARGET
#' and WEIGHT. Source and target should have GM internal ID
#' @param geneFile (char) path to GeneMANIA GENES.txt file which resolves
#' internal (GM-assigned) patient ID to names
#' @return (data.frame) network with SOURCE TARGET WEIGHT, but with 
#' original patient ID
#' @export
GM_netIDtoNames <- function(netFile,geneFile) {
		genes <- read.delim(geneFile,sep="\t",h=F,as.is=T)
		net	<- read.delim(netFile,sep="\t",h=F,as.is=T)
		
		midx <- match(net[,1],genes[,1])
		net$SOURCE <- genes[midx,2]
		midx <- match(net[,2],genes[,1])
		net$TARGET <- genes[midx,2]
		net <- net[,c(4,5,3)]
		colnames(net)[3] <- "weight"

		net
}
