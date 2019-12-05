#' Prune network by retaining strongest edges
#'
#' @param net (data.frame)  Network to prune. Columns are: source,target,weight
#' @param vertices (char) node names. Should match those in net[,seq_len(2)]
#' @param topX (numeric 0 to 1) Fraction of strongest edges to retain
#' @return (data.frame) pruned network. Three columns: AliasA, AliasB, and 
#' weight
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph delete.edges
#' @importFrom igraph get.edgelist
#' @importFrom igraph edge_attr
#' @importFrom igraph E
#' @export
pruneNetByStrongest <- function(net, vertices, topX = 0.1) {
    g <- igraph::graph_from_data_frame(net, vertices)
    weights <- sort(E(g)$weight, decreasing = TRUE)
    thresh <- weights[length(weights) * topX]
    g2 <- delete.edges(g, which(E(g)$weight < thresh))
    
    df <- as.data.frame(get.edgelist(g2))
    df[, 1] <- as.character(df[, 1])
    df[, 2] <- as.character(df[, 2])
    df$weight <- edge_attr(g2, name = "weight")
    colnames(df) <- c("AliasA", "AliasB", "weight")
    df
}

