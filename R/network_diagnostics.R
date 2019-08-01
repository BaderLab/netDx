# functions for exploring PSN

#' plot distribution of PSN edge weights
#' 
#' @details Used for comparing the properties of PSN for different similarity
#' metrics
#' @param netDir (char) path to directory of networks. Networks must be
#' in SIF format; i.e. with three columns: source, target, edge weight.
#' @param numPatients (integer) number of patients expected in network.
#' Used to calculate sparsity. If set to NULL, sparsity will not be calculated.
#' @param setName (char) used for output file name (<setName>_edgeDistr.pdf) 
#' @param randomSample (integer) if set to NULL. Plots for all networks.
#' If set to an integer plots the number of networks requested, to the max
#' number of networks available
#' @param grepFiles (char) searches for files with the provided string.
#' @param stripFromFileName (char) vector of char patterns to remove from file name
#' @param verbose (logical) print messages
#' @return No value. Side effect of creating plots of network edge-weight 
#' distribution and printing network statistics to console. Files created
#' include:
#' 1) <netDir>/<setName>_edgeDistr.pdf: Violin plots showing edge-weight
#' distribution for each network
#' 2) <netDir>/<setName>_edgeCount.pdf: Single violin plot of edge count
#' distribution across all networks.
#' The script also prints the mean and SD of edge weight across all networks
#' @importFrom caroline violins
#' @export
plotNetEdgeDistribution <- function(netDir,setName="netDir",randomSample=NULL, 
	grepFiles="cont.txt",stripFromFileName="_cont.txt",numPatients=NULL,
	verbose=TRUE) {
fList <- dir(path=netDir,pattern=grepFiles)
cat(sprintf("%i files match pattern\n", length(fList)))

if (!is.null(randomSample)) {
	randomSample <- min(randomSample, length(fList))
	cat(sprintf("randomSample set - picking %i files\n", randomSample))
	
	oldseed <- NULL
	if (exists(".Random.seed")) oldseed <- .Random.seed
	fList <- sample(fList,randomSample,FALSE)
	if (!is.null(oldseed)) .Random.seed <- oldseed
}

out <- list()
edge_count <- list()
for (curr in fList) {
    print(curr)
    dat <- read.delim(sprintf("%s/%s",netDir,curr),sep="\t",header=FALSE,as.is=TRUE)
    out[[curr]] <- dat[,3]
	edge_count[[curr]] <- nrow(dat)
}

mu <- mean(unlist(lapply(out,mean)))
med <- mean(unlist(lapply(out,median)))

pdf(sprintf("%s/%s_edgeDistr.pdf",netDir,setName),height=6,width=13)
par(mar=c(10,5,3,3),las=2)
tryCatch({
	violins(out,connect=FALSE,at=1:length(out),names=rep("",length(out)),
		main=sprintf("%s:Edge wts",setName),ylab="pairwise similarity")	
	nm <- names(out)
	for (k in 1:length(stripFromFileName)) { 
		nm <- sub(stripFromFileName[k],"",nm)
	}
	axis(1,at=1:length(out),labels=nm,cex=0.8)
},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

cat(sprintf("Network statistics\n"))
cat("----------------------------\n")
cat(sprintf("Mean edge weight = %1.2f\n",mu)) 
cat(sprintf("Median edge weight = %1.2f\n",med)) 
cat("\n")

edge_count <- unlist(edge_count)
pdf(sprintf("%s/%s_edgeCount.pdf",netDir,setName),height=9,width=5)
tryCatch({
	violins(edge_count,main=sprintf("%s: num edges",setName),
	ylab="Num. edges")
	if (!is.null(numPatients)) {
		sparsity <- ( edge_count/choose(numPatients,2) ) * 100
		violins(sparsity,main=sprintf("%s: sparsity",setName),
		ylab="% sparsity")
	}
},error=function(ex) {
	print(ex)
},finally={
	dev.off()
})

cat(sprintf("Network statistics\n"))
cat("----------------------------\n")
cat(sprintf("Edge count = %i - %i \n", min(edge_count),max(edge_count)))
cat(sprintf("Mean edge count = %1.2f\n",mean(edge_count))) 
cat(sprintf("Median edge count = %1.2f\n",median(edge_count))) 
cat("\n")
}
