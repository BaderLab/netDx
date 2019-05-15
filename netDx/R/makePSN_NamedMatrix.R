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
#' @param outDir (char) path to directory where networks are written
#' @param simMetric (char) measure of similarity. See \code{getSimilarity()}
#' for details. If writeProfiles is set to TRUE, must be one of pearson
#' (Pearson correlation) or MI (correlation by mutual information).
#' @param cutoff (numeric) patients with similarity smaller than this value
#' are not included in the corresponding interaction network
#' @param verbose (logical) print detailed messages
#' @param numCores (integer) number of cores for parallel network generation
#' @param useExtLib (logical) if TRUE, uses NetPreProc::Sparsify.matrix()
#' with user-supplied k; else uses the implementation in this file
#' @param writeProfiles (logical) use GeneMANIA's ProfileToNetworkDriver to
#' create interaction networks. If TRUE, this function writes subsets 
#' of the original data corresponding to networks to file (profiles). 
#' If FALSE, uses  getSimilarity() and writes interaction networks.
#' @param useSparsify2 (logical). Cleaner sparsification routine. 
#' If FALSE, uses new matrix-based sparsify3
#' @param sparsify_edgeMax (numeric). 
#' @param minMembers (integer) min number of measures in a network for 
#' the network to be included. Useful when similarity measures require a minimum
#' number of measures to be meaningful (e.g. minimum of 6 for Pearson correlation)
#' @param append (logical) if TRUE does not overwrite netDir.
#' @param ... passed to \code{getSimilarity()}
#' @return (char) Basename of files to which networks are written.  
#' Side effect of writing interaction networks in \code{outDir}
#' @examples data(TCGA_mini,pathwayList); 
#' # you may get a warning message that the output directory already
#' # exists; ignore it
#' out <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList, 
#' 	".",writeProfiles=TRUE)
#' @export
makePSN_NamedMatrix <- function(xpr, nm, namedSets, outDir,
	simMetric="pearson",verbose=TRUE,
	numCores=1L,writeProfiles=TRUE,
	sparsify=FALSE,useSparsify2=FALSE,cutoff=0.3,sparsify_edgeMax=1000,
	sparsify_maxInt=50,minMembers=1L,
	append=FALSE,...){
	if (!append) {
		if (file.exists(outDir)) unlink(outDir,recursive=TRUE) 
		dir.create(outDir)
	} else {
		if (!file.exists(outDir)) {
				cat("You asked for append but the directory doesn't exist. Helpfully creating it\n")
				dir.create(outDir)
		}
	}

	if ((!simMetric %in% c("pearson","MI"))  & writeProfiles==TRUE) {
	print(simMetric)
		stop("writeProfiles must only be TRUE with simMetric set to pearson or MI. For all other metrics, set writeProfiles=FALSE")
	}
	
	#if (!sparsify & useSparsify2) { stop("if useSparsify=TRUE then sparsify must also be set to TRUE\n")}

	cl	<- makeCluster(numCores,outfile=sprintf("%s/makePSN_log.txt",outDir))
	registerDoParallel(cl)

	if (simMetric=="pearson") {
		cat("simMetric set to pearson; forcing minMembers to be 5.\n")
		minMembers <- 5;
	}

	# process pathways in parallel
	outFiles <- foreach (curSet=names(namedSets)) %dopar% {
		if (verbose) cat(sprintf("%s: ", curSet))
		idx <- which(nm %in% namedSets[[curSet]])
		if (verbose) cat(sprintf("%i members\n", length(idx)))

		oFile <- NULL
 		# has sufficient connections to make network
		if (length(idx)>=minMembers) {
			if (writeProfiles) {
				outFile <- sprintf("%s/%s.profile",outDir,curSet)
				write.table(t(xpr[idx,,drop=FALSE]),file=outFile,sep="\t",
							col=F,row=T,quote=F)
			} else {
				outFile <- sprintf("%s/%s_cont.txt", outDir, curSet)
				cat(sprintf("computing sim for %s\n",curSet))
				sim 	<- getSimilarity(xpr[idx,,drop=FALSE], 
										 type=simMetric,...)
				if (is.null(sim)) {
					cat(sprintf("%s: sim is null\n",curSet))
					browser()
				}
###				if (!useSparsify2) {# prepare for internal sparsifier
###					idx <- which(upper.tri(sim,diag=F))
###					ij <- matrix_getIJ(dim(sim),idx)
###	
###					# make interaction network
###					pat_pairs <- data.frame(p1=rownames(sim)[ij[,1]], 
###										p2=colnames(sim)[ij[,2]], 
###										similarity=sim[idx])
###
###					too_weak    <- which(pat_pairs[,3] < cutoff | 
###										is.na(pat_pairs[,3]))
###					if (any(too_weak)) {
###						if (verbose) 
###							cat(sprintf("\t%i weak connections\n", 
###										length(too_weak)))
###						pat_pairs <- pat_pairs[-too_weak,]
###					}
###				} else {	 # stick to sim matrix
					pat_pairs <- sim
###				}

				if (sparsify) {
					if (useSparsify2) {
					tryCatch({
					 spmat <- sparsify2(pat_pairs,cutoff=cutoff,
							EDGE_MAX=sparsify_edgeMax,
							outFile=outFile,maxInt=sparsify_maxInt)
					},error=function(ex) {
						cat("sparse caught error\n"); browser()
					})
					} else {
						cat("sparsify3\n")
					tryCatch({
				     sp_t0 <- Sys.time()
					 spmat <- sparsify3(pat_pairs,cutoff=cutoff,
							EDGE_MAX=sparsify_edgeMax,
							outFile=outFile,maxInt=sparsify_maxInt)
					 print(Sys.time()-sp_t0)
					},error=function(ex) {
						cat("sparse caught error\n"); browser()
					})
					}
				} else {
				write.table(pat_pairs, file=outFile,sep="\t",
					col=FALSE,row=FALSE,quote=FALSE)
				print(basename(outFile))
				cat("done\n")
				}
			}
#cat("got here\n")
			oFile <- basename(outFile)
		}
		oFile
#cat("out of loop\n")
	}
	stopCluster(cl)
	outFiles
}
