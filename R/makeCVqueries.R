#' Randomly select patients for query in cross-validation
#' 
#' @param incPat (char) vector of patient IDs to be included in query
#' @param nFold (integer) n for n-fold cross validation
#' @param setSeed (integer) Set this to an integer to make results
#' 	reproducible. If NULL, seed is not explicitly set.
#' @param verbose (logical) print messages
#' @return (list) of length \code{nFold}, containing names of patients in
#' query file for each fold
#' @export
makeCVqueries <- function(incPat, nFold=10L,setSeed=42L,verbose=TRUE) {

if (!is.null(setSeed)) {
	cat(sprintf("Setting seed for reproducibility: %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}

# randomly reorder for N-fold partitioning.
incPat <- sample(incPat,replace=FALSE); 

# num in query file
num2samp	<- floor(((nFold-1)/nFold)*length(incPat))
# num to retrieve from GM database in each iteration
csize	<- round((1/nFold)*length(incPat))

if (verbose) {
	cat(sprintf("Read %i IDs\n\t%i-fold CV\n", length(incPat),nFold))
	cat(sprintf("\tEach iter will sample %i records, %i will be test\n",
				num2samp,csize))
}

out <- list()
for (k in 1:nFold) {
	sidx	<- ((k-1)*csize)+1;
	eidx	<- k*csize; 
	if (k==nFold) eidx <- length(incPat)
	if (verbose) 
			cat(sprintf("chunk %i: %i test (%i-%i); ", 
						k, eidx-sidx+1, sidx,eidx))

	out[[k]] <- setdiff(incPat, incPat[sidx:eidx])
	if (verbose) cat(sprintf("%i query\n", length(out[[k]])))
}

out
}
