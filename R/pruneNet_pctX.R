#' Prune network by retaining strongest edges
#'
#' @param net (data.frame)  Network to prune. Columns are: source,target,weight
#' @param vertices (char) node names. Should match those in net[,1:2]
#' @param pctX (numeric 0 to 1) Fraction of top/bottom edges to retain
#' @param useTop (logical) if TRUE prunes to top pctX edges; else
#' prunes to bottom pctX edges
#' @return (data.frame) pruned network. Three columns: AliasA, AliasB, and 
#' weight
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph delete.edges
#' @importFrom igraph get.edgelist
#' @importFrom igraph edge_attr
#' @importFrom igraph E
#' @export
pruneNet_pctX <- function(net,vertices, pctX=0.1, useTop=TRUE) {
	g <- igraph::graph_from_data_frame(net,vertices=vertices)
	wt <- sort(E(g)$weight, decreasing=TRUE)

	if (useTop) { # keep topmost edges
		thresh <- wt[length(wt) * pctX]
		g2 <- delete.edges(g,which(E(g)$weight < thresh))

	} else { # keep bottom-most edges 
		thresh <- wt[length(wt) * (1-pctX)]
		g2 <- delete.edges(g,which(E(g)$weight > thresh))
	}

	df <- as.data.frame(get.edgelist(g2))
	df[,1] <- as.character(df[,1])
	df[,2] <- as.character(df[,2])
	df$weight <- edge_attr(g2,name="weight")
	colnames(df) <- c("AliasA","AliasB","weight")

return(df)
}

