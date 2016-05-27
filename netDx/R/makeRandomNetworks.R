#' Generates random binary interaction networks
#' 
#' @details Used as a control to assess the robustness of GeneMANIA to 
#' random cliques.
#'
#' @param id (char) vector of patient IDs
#' @param numNets (integer) number of networks to generate
#' @param minNodes (integer) min number of patients in networks
#' @param maxNodes (integer) max number of patients in networks
#' @param sizeProb (char) Probability distribution for number of patients
#' in a networks. Options are "unif" for uniform distribution (networks
#' of all sizes between minNodes and maxNodes are equally likely) sand "exp"
#' for exponentially decaying function (more smaller networks than larger
#' networks)
#' @param netName (char) prefix for network name. 
#' @param outDir (char) output directory
#' @return char vector of network names
#' @export
makeRandomNetworks <- function(id, numNets=100L, minNodes=2L, 
	maxNodes=20L,sizeProb="exp",netName="RANDOM",outDir=".",setSeed=42L) {

if (!is.null(setSeed)) {
	cat(sprintf("Setting seed for reproducibility: %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}

# network size exponentially decays from min to max
# do this by inverse sampling
net_size <- switch(sizeProb,
	exp={
		rndsmp	<- runif(numNets)	# first sample from U[0,1]
		qx		<- qexp(rndsmp)		# map to values on exponential dist.
		a <- minNodes; b <- maxNodes
		# rescale to [minNodes, maxNodes] range
		net_size	<- ((qx/max(qx))*(b-a))+minNodes	
		net_size	<- round(net_size)
	}, unif=round(runif(numNets, minNodes, maxNodes)),
	stop("invalid value for sizeProb")
)
print(table(net_size))

outNames <- character()
for (i in 1:numNets) {
	# network i contains random sample of patients from the list
	p	<- sample(id, net_size[i],FALSE)
	pat_pairs <- t(combn(p,2))
	pat_pairs <- cbind(pat_pairs,1);

	nm	<- sprintf("%s_%i", netName,i)
	write.table(pat_pairs, 
				file=sprintf("%s/%s_cont.txt", outDir,nm),
				sep="\t",col=FALSE,row=FALSE,quote=FALSE)

	outNames <- cbind(outNames, sprintf("%s_cont.txt",nm))
}
names(outNames) <- NULL
outNames
}
