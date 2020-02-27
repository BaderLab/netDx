#' Visualize integrated patient similarity network based on selected features
#'
#' @details Generates a Cytoscape network where nodes are patients and edges
#' are weighted by aggregate pairwise patient similarity.
#' Integrated PSN plotting is intended to run
#' after feature selection, which identifies the subset of input networks
#' predictive for each class of interest. The method of generating the network
#' is as follows:
#' All networks feature-selected in either patient ategory are concatenated;
#' where a network is feature-selected in both categories, it is included once.
#' The similarity between two patients in the integrated network is the mean of
#' corresponding pairwise similarities. Dissimilarity is defined as
#' 1-similarity, and Dijkstra distances are computed on this resulting network.
#' For visualization, only edges representing the top fraction of distances
#' (strongest edge weights) are included.
#' @param pheno (data.frame) sample metadata, requires ID and STATUS column.
#' The status column will be used to retrieve patient class names
#' @param netDir (char) path to directory with networks
#' @param topX (numeric between 0 and 1) fraction of strongest edges to keep
#' e.g. topX=0.2 will keep 20\% top edges
#' @param aggFun (char) function to aggregate edges from different PSN
#' @param outDir (char) path to directory for tmp work
#' @param calcShortestPath (logical) if TRUE, computes weighted shortest path
#' Unless you plan to analyse these separately from looking at the shortest 
#' path violin plots or integrated PSN in Cytoscape, probably good to set to 
#' FALSE.
#' for all pairwise classes
#' @param verbose (logical) print detailed messages
#' @return (list) information about the integrated network
#  1) aggPSN_Full: (char) path to aggregated patient similarity network. This
#' is the unpruned network and not the one displayed in Cytoscape. Edges are
#' weighted by similarity
#' 2) aggPDN_pruned: (char) path to file with pruned patient dissimilarity
#' network, derived by inverting edge weights of the net in aggPSN_Full and
#' keeping the strongest edges.
#' 2) incNets: (char) vector of networks used for integration
#' 3) network_suid: (char) Cytoscape network identifier, so user may further
#' 4) netPng: (char) path to png file with patient dissimilarity network
#' created by Cytoscape
#' @importFrom RColorBrewer brewer.pal
#' @import RCy3
#' @export
plotIntegratedPSN <- function(setName="predictor",pheno,netDir,
	topX=0.05, aggFun="MAX",outDir=".",calcShortestPath=FALSE,
	verbose=FALSE,...) {

if (missing(pheno)) stop("pheno is missing.")
if (missing(netDir)) stop("netDir is missing.")

predClasses <- unique(pheno$STATUS)
colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
message(sprintf("%i classes: { %s }", 
	length(predClasses), paste(predClasses,
	collapse=",")))

# read patient and network identifiers
pid		<- read.delim(sprintf("%s/GENES.txt",netDir),
	,sep="\t",header=FALSE,as.is=TRUE)[,1:2]
colnames(pid)[1:2] <- c("GM_ID","ID")
netid	<- read.delim(sprintf("%s/NETWORKS.txt",netDir),
	sep="\t",header=FALSE,as.is=TRUE)
	colnames(netid)[1:2] <- c("NET_ID", "NETWORK")

# aggregate
message("* Computing aggregate net")
out <- writeWeightedNets(pid,
	netIDs=netid,netDir=sprintf("%s/INTERACTIONS",netDir),
	filterEdgeWt=0,limitToTop=Inf,
	plotEdgeDensity=FALSE,
	aggNetFunc=aggFun,verbose=FALSE)

aggNet <- out$aggNet

# calculate shortest paths among and between classes
if (calcShortestPath) {
	x <- compareShortestPath(aggNet, pheno,verbose=FALSE,plotDist=TRUE)

	gp <- unique(pheno$GROUP)
	oppName <- paste(gp[1],gp[2],sep="-")
	
	curDijk <- matrix(NA,nrow=1,ncol=6)
	colnames(curDijk) <- c(gp[1],gp[2],oppName,"overall",
				sprintf("p%s-Opp",gp[1]),sprintf("p%s-Opp",gp[2]))
	
	# mean shortest path		
	message("Shortest path averages &")
	message("p-values (one-sided WMW)")
	message("------------------------------------")
	for (k in 1:length(x$all)) {
		cur <- names(x$all)[k]
		idx <- which(colnames(curDijk) %in% cur)
		curDijk[1,idx] <- median(x$all[[k]])
		message(sprintf("\t%s: Median = %1.2f ", cur,curDijk[1,idx]))
		if (cur %in% gp) {
			tmp <- wilcox.test(x$all[[cur]],x$all[[oppName]],
				alternative="less")$p.value
			curDijk[4+k] <- tmp
			message(sprintf(" p(<Opp) = %1.2e\n", curDijk[4+k]))
		}
	}
	#write.table(curDijk,
	#	file=sprintf("%s/shortest_paths.txt",outDir),
	#	sep="\t",col.names=TRUE,row.names=TRUE,quote=FALSE)
	
}
message("got past shortest path computation")
aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity
aggNet <- aggNet[,1:3]

# create a pruned network for visualization
message("* Prune network")
colnames(aggNet) <- c("source","target","weight")
aggNet_pruned <- pruneNetByStrongest(aggNet,pheno$ID, topX=topX)
message("* Past pruning")

# plot in Cytoscape
message("* Creating network in Cytoscape")
colnames(pheno)[which(colnames(pheno)=="ID")] <- "id"
colnames(aggNet_pruned) <- c("source","target","weight")
# layout network in Cytoscape
createNetworkFromDataFrames(
	nodes=pheno, edges=aggNet_pruned,
	title=sprintf("%s_%s_top%1.2f",setName,aggFun,topX),
	collName=setName
)

styleName <- "PSNstyle"
# create Cytoscape style for PSN
pal <- suppressWarnings(brewer.pal(name="Dark2",n=length(predClasses)))
message("* Creating style")
nodeFills <- mapVisualProperty('node fill color',"GROUP",
									mapping.type="d",
                  table.column.values=predClasses,
									visual.prop.values=pal)
defaults <- list("NODE_SHAPE"="ellipse",
                    "NODE_SIZE"=30,
                    "EDGE_TRANSPARENCY"=120,
                    "EDGE_STROKE_UNSELECTED_PAINT"="#999999",
                    "NODE_TRANSPARENCY"=120)
sty <- createVisualStyle(styleName,
                            defaults,list(nodeFills))
setVisualStyle(styleName)

message("* Applying layout\n")
layoutNetwork("kamada-kawai column=weight")
setVisualStyle(styleName)

out <- list(
	patientSimNetwork_unpruned=aggNet,
	patientDistNetwork_pruned=aggNet_pruned
)

return(out)
}
