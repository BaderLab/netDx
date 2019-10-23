#' Visualize integrated patient similarity network based on selected features
#'
#' @details Generates a Cytoscape network where nodes are patients and edges
#' are weighted by aggregate pairwise patient similarity.
#' Integrated PSN plotting is intended to run
#' after feature selection, which identifies the subset of input networks
#' predictive for each class of interest. The method of generating the network
#' is as follows:
#' All networks feature-selected in either patient category are concatenated;
#' where a network is feature-selected in both categories, it is included once.
#' The similarity between two patients in the integrated network is the mean of
#' corresponding pairwise similarities. Dissimilarity is defined as
#' 1-similarity, and Dijkstra distances are computed on this resulting network.
#' For visualization, only edges representing the top fraction of distances
#' (strongest edge weights) are included.
#' @param setName (char) network name, used for output files and network name
#' in Cytoscape.
#' @param runCytoscape (logical) if TRUE will visualize PSN in Cytoscape
#' (Cytoscape must be launched for this to happen). If FALSE, will compute
#' and write integrated PSN to file but will not visualize the network. 
#' Set to FALSE if you do not have or want to run Cytoscape (e.g. running on
#' terminal or Docker container, where Cytoscape may be unavailable).
#' @param pheno (data.frame) sample metadata, requires ID and STATUS column.
#' The status column will be used to retrieve patient class names
#' @param baseDir (char) path to the results directory one level above class-
#' specific interaction networks.
#' Standard nested cross-validation: <dataDir>/rng1 so the program finds
#' <dataDir>/rng1/<classA>/tmp/INTERACTIONS/
#' Single round of cross-validation: <dataDir>, so the program finds
#' <dataDir>/<classA>/tmp/INTERACTIONS/
#' @param netNames (list) keys are classes, and each value is a vector with
#' the names of patient similarity networks to combine for the integrated
#' network. These should probably be feature-selected networks
#' @param topX (numeric between 0 and 1) fraction of strongest edges to keep
#' e.g. topX=0.2 will keep 20\% top edges
#' @param aggFun (char) function to aggregate edges from different PSN
#' @param outDir (char) path to directory for tmp work
#' @param calcShortestPath (logical) if TRUE, computes weighted shortest path
#' @param savePaths (logical) if TRUE, writes all pairwise shortest paths to file. 
#' Unless you plan to analyse these separately from looking at the shortest path
#' violin plots or integrated PSN in Cytoscape, probably good to set to FALSE.
#' for all pairwise classes
#' @param nodeSize (integer) Node size in Cytoscape-rendered PSN 
#' @param edgeTransparency (integer) Edge transparency in Cytoscape-rendered
#' @param edgeWidth (integer) Edge line width in Cytoscape-rendered PSN.
#' @param nodeTransparency (integer) Node transparency in Cytoscape-rendered
#' PSN.
#' @param edgeStroke (character) hex colour for edge stroke.
#' @param nodePalette (character) RColorBrewer palette for node colours.
#' @param imageFormat (character) file format to export PSN to image. One of
#' JPEG, PDF, PNG, or SVG (see RCy3::exportImage()).
#' @param verbose (logical) print detailed messages
#' @importFrom stats wilcox.test 
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
#' @examples
#' data(KIRC_pheno)
#' inDir <- sprintf("%s/extdata/example_output", path.package("netDx"))
#' outDir <- paste(getwd(),"plots",sep="/")
#' if (!file.exists(outDir)) dir.create(outDir)
#' featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))
#' featSelNet <- lapply(featScores, function(x) {
#' callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
#' })
#' plotIntegratedPSN(pheno=KIRC_pheno,baseDir=sprintf("%s/rng1",inDir),
#'	netNames=featSelNet,outDir=outDir,runCytoscape=FALSE)
#' @import httr
#' @import ggplot2
#' @import RColorBrewer
#' @import RCy3
#' @export
plotIntegratedPSN <- function(setName="predictor",pheno,baseDir,netNames,
	topX=0.2, aggFun="MEAN",outDir=tempdir(),
	calcShortestPath=TRUE,savePaths=FALSE,
	nodeSize=100,nodeTransparency=200,
	edgeTransparency=120,edgeStroke="#999999",edgeWidth=1,imageFormat="PNG",
	nodePalette="Dark2",verbose=FALSE,runCytoscape=TRUE) { 

if (missing(pheno)) stop("pheno is missing.")
if (missing(baseDir)) stop("baseDir is missing.")
if (missing(netNames)) stop("netNames is missing.")

predClasses <- unique(pheno$STATUS)
colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
message(sprintf("%i classes: { %s }\n", length(predClasses), paste(predClasses,
	collapse=",")))

# compile node/network locations for dataset
ptFile 	<- list()
netInfo <- list()
netDir	<- list()
for (cur in predClasses) {
	ptFile[[cur]] <- sprintf("%s/%s/tmp/GENES.txt",baseDir,cur)
	netInfo[[cur]]<- sprintf("%s/%s/tmp/NETWORKS.txt",baseDir,cur)
	netDir[[cur]]	<- sprintf("%s/%s/tmp/INTERACTIONS", baseDir,cur)
}


# check: patient IDs should be the same for both classes
print(ptFile)
if (length(ptFile)>=2) {
for (k in 2:length(ptFile)) {
	tmp <- system2(sprintf("diff %s %s", ptFile[[1]],ptFile[[k]]),stdout=TRUE)
	if (!is.null(attr(tmp,"status"))){
		message(sprintf("Gene IDs for %s doesn't match that for %s.\n", 
				predClasses[k],predClasses[1]))
		message("These should be identical. Please check.\n")
	}
}

# pool best nets for both
poolDir <- sprintf("%s/pool",outDir)
if (file.exists(poolDir)) unlink(poolDir,recursive=TRUE)
dir.create(poolDir)
	newNetIDs <- list()

# pool feature selected nets from both groups
alreadyAdded <- c()
for (gps in names(ptFile)) {
	message(sprintf("Group %s\n", gps))
	pTally <- netNames[[gps]]
	pTally <- sub("_cont|\\.profile","",pTally)

	netIDs <- read.delim(netInfo[[gps]],sep="\t",h=FALSE,as.is=TRUE)
	netIDs[,2] <- sub("_cont|\\.profile","",netIDs[,2])

	idx <- which(pTally %in% alreadyAdded)
	if (any(idx)) {
		if (verbose)
			message(sprintf("Found a net already added before, removing: {%s}\n",
				paste(pTally[idx],collapse=",")))
		pTally <- pTally[-idx]
	}

	if (length(pTally)>=1) alreadyAdded <- c(alreadyAdded,pTally)

	if (verbose) print(pTally)
	curNetIds <- matrix(NA,nrow=length(pTally),ncol=2)
	ctr <- 1

	# copy feature selected nets for this group
	for (cur in pTally) {
		idx <- which(netIDs[,2] == cur)
		if (length(idx)<1) {
			stop(sprintf("plotIntegratedPSN.R: %s: mapping of network identifier to name not found!\n",cur))
		}
		netID <- sprintf("1.%s.txt",netIDs[idx,1])
		tmp <- sprintf("%s.%s", gps,netID)
		file.copy(from=sprintf("%s/%s", netDir[[gps]], netID),
		to=sprintf("%s/1.%s",poolDir,tmp))

		curNetIds[ctr,] <- c(sub(".txt","",tmp),
		paste(gps,netIDs[idx,2],sep="."))
		##message(sprintf("\t%s -> %s\n", cur, tmp))
		ctr <- ctr+1
	}
	newNetIDs[[gps]] <- curNetIds
}

# write compiled net info file
newNetIDs <- do.call("rbind",newNetIDs)
netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
write.table(newNetIDs,file=netInfo_combinedF,sep="\t",
col=FALSE,row=FALSE,quote=FALSE)

# aggregate
message("* Computing aggregate net\n")
aggNetFile <- netDx::writeWeightedNets(ptFile[[1]],
	netInfo=netInfo_combinedF,
	poolDir,keepNets=newNetIDs[,2],outDir,
	filterEdgeWt=0,limitToTop=Inf,
	plotEdgeDensity=FALSE,
	writeAggNet=aggFun,verbose=FALSE)

# read in aggregated similarity net and convert to dissimilarity for net view
aggNet<- read.delim(aggNetFile,sep="\t",h=TRUE,as.is=TRUE)[,1:3]
colnames(aggNet) <- c("AliasA","AliasB","weight")

# calculate shortest paths among and between classes
if (calcShortestPath) {
	x <- compareShortestPath(aggNet, pheno,verbose=FALSE,plotDist=TRUE)

	gp <- unique(pheno$GROUP)
	oppName <- paste(gp[1],gp[2],sep="-")
	
	curDijk <- matrix(NA,nrow=1,ncol=6)
	colnames(curDijk) <- c(gp[1],gp[2],oppName,"overall",
				sprintf("p%s-Opp",gp[1]),sprintf("p%s-Opp",gp[2]))
	
	# mean shortest path		
	message("Shortest path averages &\n")
	message("p-values (one-sided WMW)\n")
	message("------------------------------------\n")
	for (k in 1:length(x$all)) {
		cur <- names(x$all)[k]
		idx <- which(colnames(curDijk) %in% cur)
		curDijk[1,idx] <- median(x$all[[k]])
		message(sprintf("\t%s: Median = %1.2f\n", cur,curDijk[1,idx]))
		if (cur %in% gp) {
			tmp <- wilcox.test(x$all[[cur]],x$all[[oppName]],
				alternative="less")$p.value
			curDijk[4+k] <- tmp
			message(sprintf(" p(<Opp) = %1.2e\n", curDijk[4+k]))
		}
	}
	write.table(curDijk,file=sprintf("%s/shortest_paths.txt",outDir),
		sep="\t",col=TRUE,row=TRUE,quote=FALSE)
	#print(t(curDijk))

	if (savePaths) {
		shortest_path <- x
		save(shortest_path,file=sprintf("%s/shortest_paths.rda",outDir))
	}
	
	}

aggNet[,3] <- 1-aggNet[,3]  # convert similarity to dissimilarity

}

# create a pruned network for visualization
aggNet_pruned <- netDx::pruneNetByStrongest(aggNet,pheno$ID, topX=topX)
outFile <- sprintf("%s/%s_prunedNet_top%1.2f.txt",outDir,setName,topX)
write.table(aggNet_pruned,file=outFile,sep="\t",col=TRUE,row=FALSE,
	quote=FALSE)


if (runCytoscape) {
st <- cytoscapePing()
if (is(st,"numeric")) { # error
	stop("Please launch Cytoscape and try again.")
}

message("* Creating network in Cytoscape\n")
# layout network in Cytoscape
netName <- sprintf("%s_%s_top%1.2f",setName,aggFun,topX)
colnames(pheno)[1] <- "id"
colnames(pheno)[which(colnames(pheno)=="GROUP")] <- "group"
colnames(aggNet_pruned)[1:2] <- c("source","target")
createNetworkFromDataFrames(nodes=pheno,
		edges=aggNet_pruned, title=netName, 
		collection=setName)

# apply style to network
pal <- suppressWarnings(brewer.pal(name=nodePalette,n=length(predClasses)))
	message("* Creating style\n")
	nodeLabels <- mapVisualProperty('node label','id','p')
	nodeFills <- mapVisualProperty('node fill color','group', 'd',predClasses,pal)
	defaults <- list("NODE_SHAPE"="ellipse",
			"NODE_SIZE"=nodeSize,
			"EDGE_TRANSPARENCY"=edgeTransparency,
		  "EDGE_WIDTH"=edgeWidth,
			"EDGE_STROKE_UNSELECTED_PAINT"=edgeStroke,
			"NODE_TRANSPARENCY"=nodeTransparency)
tryCatch({
	deleteVisualStyle("PSNstyle")
}, error=function(ex) {
	print(ex)
}, finally={
})
	createVisualStyle("PSNstyle",defaults,list(nodeLabels,nodeFills))
setVisualStyle("PSNstyle")
layoutNetwork("kamada-kawai column=weight")

imgFile <- sprintf("%s/outputPDN",outDir)


if (file.exists(imgFile)) unlink(imgFile)
exportImage(filename=imgFile,type=imageFormat)

out <- list(aggPSN_FULL=aggNetFile,aggPDN_pruned=outFile,
		network.suid=getNetworkSuid(),
		incNets=alreadyAdded,netView=imgFile)
} else {
	out <- list(aggPSN_FULL=aggNetFile,aggPDN_pruned=outFile,
		incNets=alreadyAdded)
}

return(out)
}
