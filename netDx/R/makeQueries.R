#' Randomly select patients for queries for feature selection
#' 
#' @param incPat (char) vector of patient IDs to be included in query
#' @param featScoreMax (integer) Number of times to run query, usually equal to the max 
#' score for features in the design (e.g. if featScoreMax=10, then this value is 10).
#' @param setSeed (integer) Set this to an integer to make results
#' 	reproducible. If NULL, seed is not explicitly set.
#' @param verbose (logical) print messages
#' @return (list) of length \code{featScoreMax}, containing names of patients in
#' query file for each fold
#' @examples 
#' data(TCGA_mini)
#' x <- makeQueries(pheno$ID)
#' @export
makeQueries <- function(incPat, featScoreMax=10L,setSeed=42L,verbose=TRUE) {
if (!is.null(setSeed)) {
	cat(sprintf("Feature scoring: set seed: %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}

# randomly reorder for N-fold partitioning.
incPat <- sample(incPat,replace=FALSE); 

# num in query file
num2samp	<- floor(((featScoreMax-1)/featScoreMax)*length(incPat))
# num to retrieve from GM database in each iteration
csize	<- round((1/featScoreMax)*length(incPat))

if (verbose) {
	cat(sprintf("Read %i IDs\n\t%i queries\n", length(incPat),featScoreMax))
	cat(sprintf("\tEach iter will sample %i records, %i will be test\n",
				num2samp,csize))
}

out <- list()
for (k in 1:featScoreMax) {
	sidx	<- ((k-1)*csize)+1;
	eidx	<- k*csize; 
	if (k==featScoreMax) eidx <- length(incPat)
	if (verbose) 
			cat(sprintf("chunk %i: %i test (%i-%i); ", 
						k, eidx-sidx+1, sidx,eidx))

	out[[k]] <- setdiff(incPat, incPat[sidx:eidx])
	if (verbose) cat(sprintf("%i query\n", length(out[[k]])))
}

out
}
