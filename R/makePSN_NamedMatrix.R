#' Create patient networks from full matrix of named measurements
#'
#' @details Creates patient similarity networks when full matrices of 
#' data are provided (e.g. gene expression, questionnaire results). To
#' generate networks from sparse data such as CNVs or indels, use 
#' \code{makePSN_RangeSets} instead.
#' The rows of the data matrix (xpr) must be named (nm); one network is 
#' create for each named set (namedSets). There are two options for the 
#' way in which networks are created, depending on the value of
#' \code{writeProfiles}. 
#' 1. writeProfiles=TRUE: GeneMANIA is used to generate interaction networks
#' and sparsify networks. This only works if the desired measure of
#' similarity is network-level Pearson correlation; an example is networks
#' at the level of pathways. In this case, the user does not explicitly 
#' specify a similarity measure and \code{simMetric} is ignored.
#' 2. writeProfiles=FALSE: GeneMANIA is not used to generate interaction
#' networks. Rather, netDx uses \code{simMetric} to create interaction
#' networks. Networks can be sparsified by excluding weak connections 
#' (cutoff). 
#' @param xpr (matrix) rows are measurements, columns are samples. Columns
#' must be named (patient ID)
#' @param nm (character) names for measurements corresponding to row order
#' of \code{xpr}. Must match the names in the named sets specified in
#' \code{nameSets}
#' @param namedSets (list) sets of names to be grouped together. keys are
#' set names, and networks will be named as these. values are character
#' vectors corresponding to groups of names (matching those in \code{nm})
#' that are input to network generation
#' @param outDir (char) path to directory where networks are written. 
#' If missing, is set to tempdir()
#' @param simMetric (char) measure of similarity. See \code{getSimilarity()}
#' for details. If writeProfiles is set to TRUE, must be one of pearson
#' (Pearson correlation) or MI (correlation by mutual information).
#' @param cutoff (numeric) patients with similarity smaller than this value
#' are not included in the corresponding interaction network
#' @param verbose (logical) print detailed messages
#' @param numCores (integer) number of cores for parallel network generation
#' @param writeProfiles (logical) use GeneMANIA's ProfileToNetworkDriver to
#' create interaction networks. If TRUE, this function writes subsets 
#' of the original data corresponding to networks to file (profiles). 
#' If FALSE, uses  getSimilarity() and writes interaction networks.
#' @param sparsify (logical). If TRUE, sparsifies patient similarity network.
#' See useSparsify2, sparsify_edgeMax and sparsify_maxInt
#' @param useSparsify2 (logical). Cleaner sparsification routine. 
#' If FALSE, uses new matrix-based sparsify3
#' @param sparsify_maxInt (numeric) Max num edges per node in sparsified network. 
#' @param sparsify_edgeMax (numeric) Max number of edges to include in the
#' final network
#' @param minMembers (integer) min number of measures in a network for 
#' the network to be included. Useful when similarity measures require a minimum
#' number of measures to be meaningful (e.g. minimum of 6 for Pearson correlation)
#' @param runSerially (logical) set to TRUE to create nets serially, rather 
#' than in parallel
#' @param ... passed to \code{getSimilarity()}
#' @return (char) Basename of files to which networks are written.  
#' Side effect of writing interaction networks in \code{outDir}
#' @import doParallel
#' @examples data(xpr,pheno,pathwayList); 
#' # you may get a warning message that the output directory already
#' # exists; ignore it
#' out <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList, 
#' 	".",writeProfiles=TRUE)
#' @export
makePSN_NamedMatrix <- function(xpr, nm, namedSets, outDir=tempdir(),
	simMetric="pearson",verbose=TRUE,
	numCores=1L,writeProfiles=TRUE,
	sparsify=FALSE,useSparsify2=FALSE,cutoff=0.3,sparsify_edgeMax=Inf,
	sparsify_maxInt=50,minMembers=1L,runSerially=FALSE,
	...) {

	if ((!simMetric %in% c("pearson","MI"))  & writeProfiles==TRUE) {
	print(simMetric)
		stop("writeProfiles must only be TRUE with simMetric set to pearson or MI. For all other metrics, set writeProfiles=FALSE")
	}
	
	cl	<- makeCluster(numCores,outfile=sprintf("%s/makePSN_log.txt",outDir))
	if (!runSerially) {
	registerDoParallel(cl)
	} else {
		message("running serially")
	}

	if (simMetric=="pearson") {
		message("Pearson similarity chosen - enforcing min. 5 patients per net.")
		minMembers <- 5;
	}

	# process pathways in parallel
	outFiles <- foreach (curSet=names(namedSets)) %dopar% {
		if (verbose) message(sprintf("%s: ", curSet))
		idx <- which(nm %in% namedSets[[curSet]])
		if (verbose) message(sprintf("%i members", length(idx)))

		oFile <- NULL
 		# has sufficient connections to make network
		if (length(idx)>=minMembers) {
			if (writeProfiles) {
				outFile <- sprintf("%s/%s.profile",outDir,curSet)
				write.table(t(xpr[idx,,drop=FALSE]),file=outFile,sep="\t",
							col.names=FALSE,row.names=TRUE,quote=FALSE)
			} else {
				outFile <- sprintf("%s/%s_cont.txt", outDir, curSet)
				message(sprintf("computing sim for %s",curSet))
				sim 	<- getSimilarity(xpr[idx,,drop=FALSE], 
										 type=simMetric,...)
				if (is.null(sim)) {
					stop(sprintf("makePSN_NamedMatrix:%s: similarity matrix is empty (NULL)\nCheck that there isn't a mistake in the input data or similarity method of choice.\n",curSet))
				}
					pat_pairs <- sim

				if (sparsify) {
					if (useSparsify2) {
					tryCatch({
					 spmat <- sparsify2(pat_pairs,cutoff=cutoff,
							EDGE_MAX=sparsify_edgeMax,
							outFile=outFile,maxInt=sparsify_maxInt)
					},error=function(ex) {
						stop("sparsify2 caught error\n"); 
					})
					} else {
						message("sparsify3")
					tryCatch({
				     sp_t0 <- Sys.time()
					 spmat <- sparsify3(pat_pairs,cutoff=cutoff,
							EDGE_MAX=sparsify_edgeMax,
							outFile=outFile,maxInt=sparsify_maxInt,verbose=FALSE)
					 print(Sys.time()-sp_t0)
					},error=function(ex) {
						stop("sparsify3 caught error\n"); 
					})
					}
				} else {
				write.table(pat_pairs, file=outFile,sep="\t",
					col.names=FALSE,row.names=FALSE,quote=FALSE)
				print(basename(outFile))
				message("done")
				}
			}
#message("got here\n")
			oFile <- basename(outFile)
		}
		oFile
#message("out of loop\n")
	}
	stopCluster(cl)
	outFiles
}
