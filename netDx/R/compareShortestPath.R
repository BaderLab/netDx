#' compare intra-cluster shortest distance to overall shortest distance of the network
#'
#' @details Uses Dijkstra's algorithm for weighted edges. Pairwise nodes with
#' infinite distances are excluded before computing average shortest path 
#' for a network. This function requires the igraph package to be installed.
#' The showNetDist option requires the gplots package to be installed.
#' @param net (data.frame) network on which to compute shortest path. 
#' SOURCE, TARGET, WEIGHTS. 
#' Column names are ignored but expects a header row. Distances will be 
#' computed based on the third column
#' @param pheno (data.frame) Node information. ID (node name) and GROUP
#' (cluster name)
#' @param showNetDist (logical) show distance heatmap (requires gplots 
#' library)
#' @param verbose (logical) print messages
#' @examples data(silh); 
#' colnames(silh$net)[3] <- "weight"
#' compareShortestPath(silh$net, silh$groups)
#' @return (list) keys are cluster names, and values are mean shortest path 
#' for subnetworks that contain only the edges where source and target both 
#' belong to the corresponding cluster. In addition, there is an "overall" 
#' entry for the mean shortest distance for the entire network.
#' @export
compareShortestPath <- function(net,pheno,showNetDist=FALSE,verbose=TRUE) {
	colnames(net) <- c("source","target","weight")

	if (verbose) {
		cat("Weight distribution:\n")
		print(summary(net[,3]))
	}

	.getAvgD <- function(mat) {
		tmp <- mat[upper.tri(mat,diag=FALSE)]
		idx <- which(is.infinite(tmp))
		if (any(idx)) tmp <- tmp[-idx]
		if (verbose) cat(sprintf("\tN=%i distances\n", length(tmp)))

		c(mean(tmp,na.rm=TRUE), sd(tmp,na.rm=TRUE))
	}
	
	g <- igraph::graph_from_data_frame(net, vertices=pheno$ID)
	d_overall <- igraph::shortest.paths(g,algorithm="dijkstra")
	cat(sprintf("Overall: %i nodes\n",length(pheno$ID)))
	tmp <- .getAvgD(d_overall)
	if (verbose)
		cat(sprintf("All-all shortest path = %2.3f (SD=%2.3f)\n",
					tmp[1],tmp[2]))
	if (showNetDist) {
			gplots::heatmap.2(t(d_overall),trace='none',scale='none',
			dendrogram='none',main="node-level shortest path")
	}

	cnames <- unique(pheno$GROUP)
	dset <- list()
	for (curr_cl in cnames) {
		cl <- pheno$ID[which(pheno$GROUP%in% curr_cl)]
		if (verbose) cat(sprintf("\n%s: N=%i nodes\n", curr_cl,length(cl)))

		#subgraph with intra-cluster connections
		g2 <- igraph::graph_from_data_frame(
			d=net[which(net[,1]%in%cl & net[,2]%in%cl),],
			vertices=cl)
		tmp <- igraph::shortest.paths(g2,algorithm="dijkstra")
		dset[[curr_cl]] <- .getAvgD(tmp)
		tmp <- dset[[curr_cl]]
		if (verbose) 
			cat(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f)\n", 
						curr_cl,curr_cl,tmp[1],tmp[2]))
	}

	# now repeat for all pairwise classes
	cpairs <- as.matrix(combinat::combn(cnames,2))
	for (k in 1:ncol(cpairs)) {
		type1 <- pheno$ID[which(pheno$GROUP %in% cpairs[1,k])]
		type2 <- pheno$ID[which(pheno$GROUP %in% cpairs[2,k])]
		idx <- which(net[,1] %in% type1 & net[,2] %in% type2)
		idx2 <- which(net[,1] %in% type2 & net[,2] %in% type1)
		idx <- c(idx,idx2)
		g <- igraph::graph_from_data_frame(d=net[idx,])
		tmp <- igraph::shortest.paths(g,algorithm="dijkstra")
		cur <- sprintf("%s-%s\n", cpairs[1,k],cpairs[2,k])
		dset[[cur]] <- .getAvgD(tmp)
		dset[[curr_cl]] <- .getAvgD(tmp)
		tmp <- dset[[curr_cl]]
		if (verbose) 
			cat(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f)\n", 
						cpairs[1,k],cpairs[2,k],tmp[1],tmp[2]))
	}

	dset[["overall"]] <- .getAvgD(d_overall)
	return(dset)

}
