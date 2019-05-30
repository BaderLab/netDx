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
#' data(xpr,pheno,cnv_GR)
#' x <- makeQueries(pheno$ID)
#' @export
makeQueries <- function(incPat, featScoreMax=10L,setSeed=42L,verbose=TRUE) {
if (!is.null(setSeed)) {
	if (verbose) cat(sprintf("\t\tSetting RNG seed to %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}

# randomly reorder for N-fold partitioning.
incPat <- sample(incPat,replace=FALSE); 

# num in query file
num2samp	<- floor(((featScoreMax-1)/featScoreMax)*length(incPat))
# num to retrieve from GM database in each iteration
csize	<- round((1/featScoreMax)*length(incPat))

if (verbose) {
	cat(sprintf("\t\t%i IDs; %i queries (%i sampled, %i test)\n",
		 length(incPat),featScoreMax,num2samp,csize))
}

out <- list()
for (k in 1:featScoreMax) {
	sidx	<- ((k-1)*csize)+1;
	eidx	<- k*csize; 
	if (k==featScoreMax) eidx <- length(incPat)
	if (verbose) cat(sprintf("\t\tQ%i: %i test; ",k, eidx-sidx+1))

	out[[k]] <- setdiff(incPat, incPat[sidx:eidx])
	if (verbose) cat(sprintf("%i query\n", length(out[[k]])))
}

out
}
