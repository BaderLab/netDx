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
#' @param dataList (MultiAssayExperiment) patient data & labels used as input
#' @param groupList (list) features to use to create integrated patient network.
#' Identical in structure to groupList in buildPredictor() method. This is a 
#' list of lists, where the outer list corresponds to assay (e.g. mRNA,
#' clinical) and inner list to features to generate from that datatype.
#' @param makeNetFunc (function) function to create features
#' @param setName (char) name to assign the network in Cytoscape
#' @param numCores (integer) number of cores for parallel processing
#' @param prune_pctX (numeric between 0 and 1) fraction of most/least 
#' edges to keep when pruning the integrated PSN for visualization.
#' Must be used in conjunction with useTop=TRUE/FALSE
#' e.g. Setting pctX=0.2 and useTop=TRUE will keep 20\% top edges
#' @param prune_useTop (logical) when pruning integrated PSN for visualization,
#' determines whether to keep strongest edges (useTop=TRUE) or weakest edges
#' (useTop=FALSE)
#' @param aggFun (char) function to aggregate edges from different PSN
#' @param outDir (char) path to directory for intermediate files. Useful for
#' debugging.
#' @param calcShortestPath (logical) if TRUE, computes weighted shortest path
#' Unless you plan to analyse these separately from looking at the shortest 
#' path violin plots or integrated PSN in Cytoscape, probably good to set to 
#' FALSE.
#' @param showStats (logical) if FALSE, suppresses shortest path-related 
#' stats, such as one-sided WMW test for testing shorter intra-class distances
#' @param nodeSize (integer) size of nodes in Cytoscape
#' @param edgeTransparency (integer) Edge transparency. 
#' Value between 0 and 255, with higher numbers leading to more opacity.
#' @param nodeTransparency (integer) Node transparency.
#' Value between 0 and 255, with higher numbers leading to more opacity.
#' @param verbose (logical) print detailed messages
#' @param plotCytoscape (logical) If TRUE, plots network in Cytoscape.
#' Requires Cytoscape software to be installed and running on the computer
#' when the function call is being made.
#' @return (list) information about the integrated network
#  1) patientSimNetwork_unpruned (matrix) full integrated 
#' similarity network
#' 2) patientDistNetwork_pruned (matrix) the network plotted in
#' Cytoscape. Also note that this is a dissimilarity network, 
#' so that more similar nodes have smaller edge weights
#' 3) colLegend (data.frame): legend for the patient network
#' plotted in Cytoscape. Columns are node labels (STATUS) and
#' colours (colour)
#' 6) outDir (char) value of outDir parameter
#' @importFrom RColorBrewer brewer.pal
#' @importFrom stats wilcox.test qexp density
#' @export
plotIntegratedPatientNetwork <- function(dataList,groupList,makeNetFunc,
	setName="predictor",prune_pctX=0.05, prune_useTop=TRUE,
	 aggFun="MAX",calcShortestPath=FALSE,
	showStats=FALSE,
	outDir=tempdir(),numCores=1L,nodeSize=50L,edgeTransparency=40L,
	nodeTransparency=155L,plotCytoscape=FALSE,
	verbose=FALSE) {

if (missing(dataList)) stop("dataList is missing.")

dat <- dataList2List(dataList, groupList)
pheno <- dat$pheno[,c("ID","STATUS")]

if (!file.exists(outDir)) dir.create(outDir)
profDir <- paste(outDir,"profiles",sep=getFileSep())
if (!file.exists(profDir)) dir.create(profDir)

# create input networks
pheno_id <- setupFeatureDB(pheno,outDir)

createPSN_MultiData(dataList=dat$assays,groupList=groupList,
			pheno=pheno_id,
			netDir=outDir,customFunc=makeNetFunc,numCores=numCores,
			verbose=FALSE)
convertProfileToNetworks(
		netDir=profDir,
		outDir=paste(outDir,"INTERACTIONS",sep=getFileSep())
)

predClasses <- unique(pheno$STATUS)
colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"

# read patient and network identifiers
pid		<- read.delim(paste(outDir,"GENES.txt",sep=getFileSep()),
	sep="\t",header=FALSE,as.is=TRUE)[,1:2]
colnames(pid)[1:2] <- c("GM_ID","ID")
netid	<- read.delim(paste(outDir,"NETWORKS.txt",sep=getFileSep()),
	sep="\t",header=FALSE,as.is=TRUE)
	colnames(netid)[1:2] <- c("NET_ID", "NETWORK")

# aggregate
message("* Computing aggregate net")
out <- writeWeightedNets(pid,
	netIDs=netid,netDir=paste(outDir,"INTERACTIONS",sep=getFileSep()),
	filterEdgeWt=0,limitToTop=Inf,
	plotEdgeDensity=FALSE,
	aggNetFunc=aggFun,verbose=FALSE)

aggNet <- out$aggNet
aggNet <- aggNet[,1:3]

distNet <- aggNet
distNet[,3] <- 1-distNet[,3] 

# calculate shortest paths among and between classes
if (calcShortestPath) {
	x <- compareShortestPath(distNet, 
		pheno,verbose=showStats,plotDist=TRUE)

	gp <- unique(pheno$GROUP)
	oppName <- paste(gp[1],gp[2],sep="-")
	
	curDijk <- matrix(NA,nrow=1,ncol=6)
	colnames(curDijk) <- c(gp[1],gp[2],oppName,"overall",
				sprintf("p%s-Opp",gp[1]),sprintf("p%s-Opp",gp[2]))
	
	# mean shortest path		
	if (showStats) {
		message("Shortest path averages &")
		message("p-values (one-sided WMW)")
		message("------------------------------------")
	}

	for (k in 1:length(x$all)) {
		cur <- names(x$all)[k]
		idx <- which(colnames(curDijk) %in% cur)
		curDijk[1,idx] <- median(x$all[[k]])
		if (showStats) message(sprintf("\t%s: Median = %1.2f ", cur,curDijk[1,idx]))
		if (cur %in% gp) {
			tmp <- wilcox.test(x$all[[cur]],x$all[[oppName]],
				alternative="less")$p.value
			curDijk[4+k] <- tmp
			if (showStats) message(sprintf(" p(<Opp) = %1.2e\n", curDijk[4+k]))
		}
	}
}
message("")

# create a pruned network for visualization
message("* Prune network")
colnames(aggNet) <- c("source","target","weight")
aggNet_pruned <- pruneNet_pctX(
	aggNet,pheno$ID, pctX=prune_pctX,useTop=prune_useTop
)
#aggNet_pruned <- distNet_pruned
#aggNet_pruned[,3] <- 1-distNet_pruned[,3]

# plot in Cytoscape
if (plotCytoscape) {
	if (!requireNamespace("RCy3",quietly=TRUE)) {
		stop("Package \"RCy3\" needed for this function to work. Please install it.",
		call.=FALSE)
	}
	message("* Creating network in Cytoscape")
	colnames(pheno)[which(colnames(pheno)=="ID")] <- "id"
	colnames(aggNet_pruned) <- c("source","target","weight")
	# layout network in Cytoscape
	if (prune_useTop) { topTxt <- "top" } else { topTxt <- "bottom" }
	RCy3::createNetworkFromDataFrames(
		nodes=pheno, edges=aggNet_pruned,
		title=sprintf("%s_%s_%s%1.2f",setName,aggFun,topTxt,prune_pctX),
		collName=setName
	)
	
	styleName <- "PSNstyle"
	# create Cytoscape style for PSN
	pal <- suppressWarnings(
		brewer.pal(name="Dark2",n=length(predClasses)))
	message("* Creating style")
	colLegend <- data.frame(STATUS=predClasses,colour=pal[1:length(predClasses)],
		stringsAsFactors=FALSE)
	nodeFills <- RCy3::mapVisualProperty('node fill color',"GROUP",
										mapping.type="d",
	                  table.column.values=predClasses,
										visual.prop.values=pal)
	defaults <- list("NODE_SHAPE"="ellipse",
	                    "NODE_SIZE"=nodeSize,
	                    "EDGE_TRANSPARENCY"=edgeTransparency,
	                    "EDGE_STROKE_UNSELECTED_PAINT"="#999999",
	                    "NODE_TRANSPARENCY"=nodeTransparency)
	sty <- RCy3::createVisualStyle(styleName,
	                            defaults,list(nodeFills))
	RCy3::setVisualStyle(styleName)
	
	message("* Applying layout\n")
	RCy3::layoutNetwork("kamada-kawai column=weight")
	RCy3::setVisualStyle(styleName)
	out <- list(
		patientSimNetwork_unpruned=aggNet,
		patientSimNetwork_pruned=aggNet_pruned,
		colLegend=colLegend,
		outDir=outDir
	)
} else {
	message(paste("* plotCytoscape is set to FALSE.",
		"Set to TRUE to visualize patient network in Cytoscape",sep=""))
	out <- list(
		patientSimNetwork_unpruned=aggNet,
		patientDistNetwork_pruned=aggNet_pruned,
		aggFun=aggFun,
		outDir=outDir
	)
}

return(out)
}
