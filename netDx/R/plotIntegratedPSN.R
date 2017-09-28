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
#' @examples
#' phenoFile <- sprintf("%s/extdata/KIRC_pheno.rda",path.package("netDx.examples"))
#' load(phenoFile)
#' inDir <- sprintf("%s/extdata/KIRC_output", path.package("netDx.examples"))
#' outDir <- paste(getwd(),"plots",sep="/")
#' if (!file.exists(outDir)) dir.create(outDir)
#' featScores <- getFeatureScores(inDir,predClasses=c("SURVIVEYES","SURVIVENO"))
#' featSelNet <- lapply(featScores, function(x) {
#' callFeatSel(x, fsCutoff=10, fsPctPass=0.7)
#' })
#' plotIntegratedPSN(pheno=pheno,baseDir=sprintf("%s/rng1",inDir),
#'	netNames=featSelNet,outDir=outDir)
#' @import httr
#' @import ggplot2
#' @import RColorBrewer
#' @export
plotIntegratedPSN <- function(setName="predictor",pheno,baseDir,netNames,
	topX=0.2, aggFun="MEAN",outDir=".",calcShortestPath=TRUE,savePaths=FALSE,
	verbose=FALSE,...) {

if (missing(pheno)) stop("pheno is missing.")
if (missing(baseDir)) stop("baseDir is missing.")
if (missing(netNames)) stop("netNames is missing.")

predClasses <- unique(pheno$STATUS)
colnames(pheno)[which(colnames(pheno)=="STATUS")] <- "GROUP"
cat(sprintf("%i classes: { %s }\n", length(predClasses), paste(predClasses,
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

# ------------------------
# Set up Cytoscape
styleName <- "PSNstyle"
portNum <- 1234
base.url <- sprintf("http://localhost:%i/v1",portNum)
res <- NULL
tryCatch({
	res	<- httr::GET(sprintf("%s/styles",base.url))
}, error=function(ex) {
	 # if Cytoscape isn't launched we may want to launch it.
	launchCytoscape()
}, finally={
	res	<- httr::GET(sprintf("%s/styles",base.url))
})
curStyles <- gsub("\\\"","",rawToChar(res$content))
curStyles <- unlist(strsplit(curStyles,","))

# throws warning if n < 3, ignore
pal <- suppressWarnings(brewer.pal(name="Dark2",n=length(predClasses)))

# create Cytoscape style for PSN
if (any(grep("PSNstyle",curStyles))) {
	cat("* Style exists, not creating\n")
} else {
	cat("* Creating style\n")
	nodeFills <- EasycyRest::map_NodeFillDiscrete("GROUP",predClasses,pal)
	defaults <- list("NODE_SHAPE"="ellipse",
			"NODE_SIZE"=30,
			"EDGE_TRANSPARENCY"=120,
			"EDGE_STROKE_UNSELECTED_PAINT"="#999999",
			"NODE_TRANSPARENCY"=120)
	sty <- createStyle(styleName,
		defaults=defaults,
		mappings=list(nodeFills))
}

# check: patient IDs should be the same for both classes
if (length(ptFile)>=2) {
for (k in 2:length(ptFile)) {
	tmp <- system(sprintf("diff %s %s", ptFile[[1]],ptFile[[k]]),intern=TRUE)
	if (!is.null(attr(tmp,"status"))){
		cat(sprintf("Gene IDs for %s doesn't match that for %s.\n", 
				predClasses[k],predClasses[1]))
		cat("These should be identical. Please check.\n")
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
	cat(sprintf("Group %s\n", gps))
	pTally <- netNames[[gps]]
	pTally <- sub("_cont|\\.profile","",pTally)

	netIDs <- read.delim(netInfo[[gps]],sep="\t",h=F,as.is=T)
	netIDs[,2] <- sub("_cont|\\.profile","",netIDs[,2])

	idx <- which(pTally %in% alreadyAdded)
	if (any(idx)) {
		if (verbose)
			cat(sprintf("Found a net already added before, removing: {%s}\n",
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
			cat(sprintf("%s: index not found!\n",cur))
			browser()
		}
		netID <- sprintf("1.%s.txt",netIDs[idx,1])
		tmp <- sprintf("%s.%s", gps,netID)
		file.copy(from=sprintf("%s/%s", netDir[[gps]], netID),
		to=sprintf("%s/1.%s",poolDir,tmp))

		curNetIds[ctr,] <- c(sub(".txt","",tmp),
		paste(gps,netIDs[idx,2],sep="."))
		##cat(sprintf("\t%s -> %s\n", cur, tmp))
		ctr <- ctr+1
	}
	newNetIDs[[gps]] <- curNetIds
}

# write compiled net info file
newNetIDs <- do.call("rbind",newNetIDs)
netInfo_combinedF <- sprintf("%s/netInfo.txt", poolDir)
write.table(newNetIDs,file=netInfo_combinedF,sep="\t",
col=F,row=F,quote=F)

# aggregate
cat("* Computing aggregate net\n")
aggNetFile <- netDx::writeWeightedNets(ptFile[[1]],
	netInfo=netInfo_combinedF,
	poolDir,keepNets=newNetIDs[,2],outDir,
	filterEdgeWt=0,limitToTop=Inf,
	plotEdgeDensity=FALSE,
	writeAggNet=aggFun,verbose=FALSE)

# read in aggregated similarity net and convert to dissimilarity for net view
aggNet<- read.delim(aggNetFile,sep="\t",h=T,as.is=T)[,1:3]
colnames(aggNet) <- c("AliasA","AliasB","weight")

# calculate shortest paths among and between classes
if (calcShortestPath) {
	pdf(sprintf("%s/shortest_paths.pdf",outDir),width=11,height=6)
	tryCatch({
		x <- compareShortestPath(aggNet, pheno,verbose=FALSE,plotDist=TRUE)
	},error=function(ex){
		print(ex)
	}, finally={
		dev.off()
	})

	gp <- unique(pheno$GROUP)
	oppName <- paste(gp[1],gp[2],sep="-")
	
	curDijk <- matrix(NA,nrow=1,ncol=6)
	colnames(curDijk) <- c(gp[1],gp[2],oppName,"overall",
				sprintf("p%s-Opp",gp[1]),sprintf("p%s-Opp",gp[2]))
	
	# mean shortest path		
	cat("Shortest path averages &\n")
	cat("p-values (one-sided WMW)\n")
	cat("------------------------------------\n")
	for (k in 1:length(x$all)) {
		cur <- names(x$all)[k]
		idx <- which(colnames(curDijk) %in% cur)
		curDijk[1,idx] <- median(x$all[[k]])
		cat(sprintf("\t%s: Median = %1.2f ", cur,curDijk[1,idx]))
		if (cur %in% gp) {
			tmp <- wilcox.test(x$all[[cur]],x$all[[oppName]],
				alternative="less")$p.value
			curDijk[4+k] <- tmp
			cat(sprintf(" p(<Opp) = %1.2e\n", curDijk[4+k]))
		}
	}
	write.table(curDijk,file=sprintf("%s/shortest_paths.txt",outDir),
		sep="\t",col=T,row=T,quote=F)
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

cat("* Creating network in Cytoscape\n")
# layout network in Cytoscape
network.suid <- EasycyRest::createNetwork(
	nodes=pheno, nodeID_column="ID",edges=aggNet_pruned,
	netName=sprintf("%s_%s_top%1.2f",setName,aggFun,topX),
	collName=setName
)
cat("* Applying layout\n")
# spring-embedded layout on edge 'weight' column
minwt <- min(aggNet_pruned$weight)
maxwt <- max(aggNet_pruned$weight)
layout.url <- sprintf("%s/apply/layouts/kamada-kawai/%s?column=weight",
	base.url,network.suid, sep="/")
response <- httr::GET(url=layout.url)
# apply style
cat("* Applying style\n")
apply.style.url <- sprintf("%s/apply/styles/%s/%i",
	base.url,styleName,network.suid)
response <- httr::GET(apply.style.url)

# fit content
fitCommand <- sprintf("%s/commands/view/fit content",base.url)
response	<- httr::GET(fitCommand)

# export to png
cat("* Exporting to PNG\n")
pngFile 		<- sprintf("%s/outputPDN.png",outDir)
if (file.exists(pngFile)) unlink(pngFile) # avoid the "overwrite file?"
																						# dialog
	exportURL <- sprintf("%s/commands/view/export?OutputFile=%s",
			base.url,pngFile)
#print(exportURL)
	response	<- httr::GET(exportURL)
#exportImage(pngFile,"PNG")

out <- list(aggPSN_FULL=aggNetFile,aggPDN_pruned=outFile,
		incNets=alreadyAdded,network_suid=network.suid,
		netView=pngFile)

return(out)
}
