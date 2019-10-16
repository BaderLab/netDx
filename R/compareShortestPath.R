#' compare intra-cluster shortest distance to overall shortest distance of the
#'  network
#'
#' @details Uses Dijkstra's algorithm for weighted edges. Pairwise nodes with
#' infinite distances are excluded before computing average shortest path 
#' for a network. This function requires the igraph package to be installed.
#' @param net (data.frame) network on which to compute shortest path. 
#' SOURCE, TARGET, WEIGHTS. 
#' Column names are ignored but expects a header row. Distances will be 
#' computed based on the third column
#' @param pheno (data.frame) Node information. ID (node name) and GROUP
#' (cluster name)
#' @param plotDist (logical) if TRUE, creates a violin plot showing the 
#' shortest path distributions for each group.
#' @param verbose (logical) print messages
#' @import ggplot2
#' @importFrom combinat combn
#' @examples data(silh); 
#' colnames(silh$net)[3] <- "weight"
#' compareShortestPath(silh$net, silh$groups)
#' @return (list) Two lists, "avg" and "all". keys are cluster names. 
#' values for "avg" are mean shortest path ; for "all", are all pairwise
#' shortest paths
#' for subnetworks that contain only the edges where source and target both 
#' belong to the corresponding cluster. In addition, there is an "overall" 
#' entry for the mean shortest distance for the entire network.
#' @importFrom stats sd
#' @export
compareShortestPath <- function(net,pheno, plotDist=FALSE,
	verbose=TRUE){	
	colnames(net) <- c("source","target","weight")

	if (verbose) {
		message("Weight distribution:\n")
		print(summary(net[,3]))
	}

	.getAvgD <- function(mat) {
		tmp <- mat[upper.tri(mat,diag=FALSE)]
		idx <- which(is.infinite(tmp))
		if (any(idx)) tmp <- tmp[-idx]
		##if (verbose) message(sprintf("\tN=%i distances\n", length(tmp)))

		c(mean(tmp,na.rm=TRUE), sd(tmp,na.rm=TRUE),length(tmp))
	}

	.getAllD <- function(mat) {
		tmp <- mat[upper.tri(mat,diag=FALSE)]
		idx <- which(is.infinite(tmp))
		if (any(idx)) tmp <- tmp[-idx]
		tmp
	}
	
	g <- igraph::graph_from_data_frame(net, vertices=pheno$ID)
	d_overall <- igraph::shortest.paths(g,algorithm="dijkstra")
	message(sprintf("Overall: %i nodes\n",length(pheno$ID)))
	tmp <- .getAvgD(d_overall)
	if (verbose)
		message(sprintf("All-all shortest path = %2.3f (SD=%2.3f) (N=%i distances)\n",
					tmp[1],tmp[2],tmp[3]))

	cnames <- unique(pheno$GROUP)
	dset <- list()
	dall <- list()
	for (curr_cl in cnames) {
		cl <- pheno$ID[which(pheno$GROUP%in% curr_cl)]
		if (verbose) message(sprintf("\n%s: N=%i nodes\n", curr_cl,length(cl)))

		#subgraph with intra-cluster connections
		g2 <- igraph::graph_from_data_frame(
			d=net[which(net[,1]%in%cl & net[,2]%in%cl),],
			vertices=cl)
		tmp <- igraph::shortest.paths(g2,algorithm="dijkstra")
		message(sprintf("%s\n",curr_cl))
		dset[[curr_cl]] <- .getAvgD(tmp)
		dall[[curr_cl]] <- .getAllD(tmp)
		tmp <- dset[[curr_cl]]
		if (verbose) 
			message(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f) (N=%i dist)\n", 
						curr_cl,curr_cl,tmp[1],tmp[2],tmp[3]))
	}

	# now repeat for all pairwise classes
	cpairs <- as.matrix(combinat::combn(cnames,2))
	message("Pairwise classes:\n")
	for (k in 1:ncol(cpairs)) {
		type1 <- pheno$ID[which(pheno$GROUP %in% cpairs[1,k])]
		type2 <- pheno$ID[which(pheno$GROUP %in% cpairs[2,k])]
		idx <- which(net[,1] %in% type1 & net[,2] %in% type2)
		idx2 <- which(net[,1] %in% type2 & net[,2] %in% type1)
		idx <- c(idx,idx2)
		g <- igraph::graph_from_data_frame(d=net[idx,])
		tmp <- igraph::shortest.paths(g,algorithm="dijkstra")
		cur <- sprintf("%s-%s", cpairs[1,k],cpairs[2,k])
		dset[[cur]] <- .getAvgD(tmp)
	#	dset[[curr_cl]] <- .getAvgD(tmp)
		dall[[cur]] <- .getAllD(tmp)
		
		tmp <- dset[[curr_cl]]
		if (verbose) 
			message(sprintf("\t%s-%s: Mean shortest = %2.3f (SD= %2.3f) (N=%i dist)\n", 
						cpairs[1,k],cpairs[2,k],tmp[1],tmp[2],tmp[3]))
	}

	dset[["overall"]] <- .getAvgD(d_overall)
	dall[["overall"]] <- .getAllD(d_overall)

	out <- list(avg=dset,all=dall)
	if (plotDist) {	
		par(las=1,bty='n')
		dl <- data.frame(intType=rep(names(dall),lapply(dall,length)),
				dijk=unlist(dall))
		plotList <- list()
		p <-ggplot(dl,aes(intType, dijk))
		p <- p + ylab("Pairwise Dijkstra distance\n(smaller is better)") 
		p <- p + xlab("Pair groups")
		p2 <- p+geom_violin(scale="width")+geom_boxplot(width=0.02) 
		print(p2)
		out[["plot"]] <- p2
		}
	return(out)

}
