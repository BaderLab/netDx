#' Convert a network in source-target-weight format to symmetric matrix
#'
#' @details A common format for network representation is to use a three
#' column table listing source node, target node, and weight.  
#' This is the format netDx uses for network integration and visualization
#' in Cytoscape. However, some functionality requires a square symmetric
#' adjacency matrix. This function takes as input the three-column format
#' and converts to the adjacency matrix. 
#' NOTE: Symmetric attribute is assumed, and the function automatically sets
#' a[i,j] = a[j,i]. Diagonal is assumed to have value of 1.0. Finally
#' missing edges will be assigned NA values.
#' @param x (data.frame) three columns, with source node, target node, and 
#' edge weight. Entries must include universe of nodes; those with missing
#' edges must be included as having edge weight NA
#' @param verbose (logical) print messages
#' @return (matrix) symmetric adjacency matrix
#' @examples
#' src <- c("A","B"); tgt <- c("C","C")
#' cur <- data.frame(source=src,target=tgt,weight=c(0.3,0.8))
#' makeSymmetric(cur)
#' @export
makeSymmetric <- function(x,verbose=FALSE) {
samps <- unique(c(x[,1],x[,2]))
newmat <- matrix(NA, nrow=length(samps),ncol=length(samps))
rownames(newmat) <- samps
colnames(newmat) <- samps
i <- 1
for (k in samps) {
	idx <- which(x[,1] == k)
	if (verbose) message(k)
	for (curr in idx) {
		#message(paste("\t",x[curr,2]))
		j <- which(colnames(newmat) == x[curr,2])
		newmat[i,j] <- x[curr,3]
		newmat[j,i] <- x[curr,3]
	}
	i <- i+1
}

diag(newmat) <- 1
return(newmat)
}

