
#' shuffle gene-to-pathway mapping
#' @param (list) keys are pathways, values are genes
#' @param (integer) rng seed
shufflePathways <- function(pathList,setSeed=42) {
cat("shuffling gene-pathway mapping\n")
t0 <- Sys.time()
set.seed(setSeed)
# count genes in pathway
ln		<- unlist(lapply(pathList, length))
genes	<- sample(unlist(pathList),replace=F) # shuffle
newP <- list()
ctr <- 1
for (k in 1:length(ln)) {
	epos <- ctr+(ln[k]-1)
	newP[[names(pathList)[k]]] <- genes[ctr:epos]
	ctr <- ctr+ln[k]
}
t1 <- Sys.time()
print(t1-t0)
return(newP)
}
