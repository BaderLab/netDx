#' Randomly select patients for queries for feature selection
#' 
#' @param incPat (char) vector of patient IDs to be included in query
#' @param featScoreMax (integer) Number of times to run query, usually equal to the max 
#' score for features in the design (e.g. if featScoreMax=10, then this value is 10).
#' @param verbose (logical) print messages
#' @return (list) of length \code{featScoreMax}, containing names of patients in
#' query file for each fold
#' @examples 
#' data(xpr,pheno,cnv_GR)
#' x <- makeQueries(pheno$ID)
#' @export
makeQueries <- function(incPat, featScoreMax=10L,verbose=TRUE) {

# randomly reorder for N-fold partitioning.
incPat <- sample(incPat,replace=FALSE); 
# num in query file
num2samp	<- floor(((featScoreMax-1)/featScoreMax)*length(incPat))
# num to retrieve from GM database in each iteration
csize	<- round((1/featScoreMax)*length(incPat))

if (verbose) {
	message(sprintf("\t\t%i IDs; %i queries (%i sampled, %i test)",
		 length(incPat),featScoreMax,num2samp,csize))
}

out <- list()
for (k in 1:featScoreMax) {
	sidx	<- ((k-1)*csize)+1;
	eidx	<- k*csize; 
	if (k==featScoreMax) eidx <- length(incPat)
	p1 <- sprintf("\t\tQ%i: %i test; ",k, eidx-sidx+1)

	out[[k]] <- setdiff(incPat, incPat[sidx:eidx])
	if (verbose) message(sprintf("%s %i query", p1, length(out[[k]])))
}

out
}
