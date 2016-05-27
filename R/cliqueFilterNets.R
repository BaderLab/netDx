#' Score networks based on their edge bias towards (+,+) interactions
#' 
#' @details Determines which networks are statistically enriched for 
#' interactions between the class of interest. The resulting \code{ENR} 
#' score and corresponding p-value serve as a filter to exclude random-like
#' interaction networks before using feature selection. This filter is
#' known to be important when patient networks are sparse and binary; e.g.
#' networks based on shared overlap of CNV locations.  If the filter is 
#' not applied, GeneMANIA WILL promote networks with slight bias towards 
#' (+,+) edges , even if these are small and random-like.
#' 
#' The measure of (+,+)-enrichment is defined as: 
#' ENR(network N) = ((num (+,+) edges) - (num other edges))/(num edges).
#' A p-value for per-network ENR is obtained non-parametrically by
#' measuring a null distribution for ENR following multiple permutations
#' of case-control labels.
#' @param netDir (char) path to dir containing all networks
#' @param pheno_DF (data.frame) for details see \code{getEnr()}
#' @param outDir (char) path to dir where output/log files are written 
#' @param numReps (integer) Max num reps for shuffling class status. Adaptive permutation is
#' used so in practice, few networks would be evaluated to this extent
#' @param minEnr (numeric from -1 to 1) Only include networks with ENR
#'	value greater than this threshold. 
#' @param outPref (char) prefix for log file (not counting the dir name)
#' @param verbose (logical) print messages
#' @param setSeed (integer) if not NULL, integer is set as seed
#' 	to ensure reproducibility in random number generation
#' @param enrType (char) see getEnr()
#' @param numCores (integer) num cores for parallel ENR computation of
#' all networks
#' @param predClass (char) see \code{getEnr()}
#' @param tmpDir (char) path to dir where temporary work can be stored
#' @param ... parameters for \code{countIntType_batch()}. 
#' @return (data.frame) networks stats from clique-filtering, one record per network
#' @export
cliqueFilterNets <- function(netDir,pheno_DF,outDir,numReps=50L,
	minEnr=-1,outPref="cliqueFilterNets",verbose=TRUE,setSeed=42L,
	enrType="binary",numCores=1L,predClass,tmpDir="/tmp",...) {
		
today	<- format(Sys.Date(),"%y%m%d")
runtime	<- format(Sys.time(),"%H%M")

# ---------------------------------------------------------------------
# work begins
logFile	<- sprintf("%s/%s_%s.log", outDir,outPref,today)
sink(logFile,split=TRUE)

tryCatch({

cat(sprintf("Network dir:\n\t%s\n", netDir))
cat(sprintf("Output dir:\n\t%s\n", outDir))
cat("\n")

if (!is.null(setSeed)) {
	cat(sprintf("Setting seed for reproducibility: %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}
cl	<- makeCluster(numCores)
registerDoParallel(cl)

# get enrichment for real networks
orig_enr <- getEnr(netDir,pheno_DF,predClass=predClass,tmpDir=tmpDir,
				   enrType=enrType,...)
plusID		<- orig_enr[["plusID"]]
minusID		<- orig_enr[["minusID"]]
orig		<- orig_enr[["orig"]]
orig_rat	<- orig_enr[["orig_rat"]]
fList		<- orig_enr[["fList"]]

print(summary(orig_rat))

idx		<- which(orig_rat >= minEnr)
cat(sprintf("\n\t%i of %i networks have ENR >= %1.1f -> filter\n", 
	length(idx), length(fList), minEnr))
fList		<- fList[idx]
orig_rat	<- orig_rat[idx]
orig		<- orig[idx,]

# now shuffle
cat("* Computing shuffled\n")
both	<- c(plusID,minusID)
n1		<- length(plusID); n <- length(both)

# index of networks that are clearly going to be not significant
# and for which it is useless to permute further 
drop_out	<- integer()
N			<- length(fList)
to_run		<- 1:N
currRep		<- 1

# rows are networks, columns are (pp-mp)/(pp+mp) for each shuffled rep
shuf_rat	<- matrix(NA,nrow=length(fList),ncol=numReps)
cat("\t")

x0 <- system.time( 
while ((length(to_run)>0) &  (currRep <=numReps)) {
	shuf	<- sample(both,replace=F) # shuffle case-control label
	# recompute pp and pm for each network
	# this step is run in parallel
	tmp		<- countIntType_batch(fList[to_run],
				shuf[1:n1], shuf[(n1+1):n],tmpDir=tmpDir,
				enrType=enrType,...)
	if (length(to_run)<2) tmp <- matrix(tmp,ncol=2)
		
	if (enrType=="binary"){
		shuf_rat[to_run,currRep]<- (tmp[,1]-tmp[,2])/(tmp[,1]+tmp[,2])
	} else if (enrType=="corr") {
		# divide by two to rescale to [-1,1]
		shuf_rat[to_run,currRep]<- (tmp[,1]-tmp[,2])/2
	} else {
		shuf_rat[to_run,currRep]<- NA
	}

	# every 10th run, identify networks where the shuffle does
	# as well as the real data >=50% of the time;
	# add these to drop-out and don't evaluate them in the future
	if (currRep %%10 == 0) {
		cat(".")
		orig_pct	<-  sapply(to_run, function(k) {
			return(sum(shuf_rat[k,1:currRep]>= orig_rat[k])/currRep)
		})

		y	<- exp(-(currRep-1)/100) # decay term
		y2	<- 3/currRep
		if (y - y2 > 0) chk_thresh <- y else chk_thresh <- y2
		idx	<- which(orig_pct >= chk_thresh)
		if (length(idx)>0)	
			drop_out	<- c(drop_out, to_run[idx])
		to_run		<- setdiff(1:N,drop_out)

		if (verbose) cat(sprintf("%i (%1.5f) %i drop out; %i left\n", 
			currRep, chk_thresh, length(idx), length(to_run)))
	}
	if (currRep%%100==0 & !verbose) cat("\n\t")
	currRep		<- currRep + 1
}
)

stopCluster(cl)

# consolidate results
cat("* Consolidating results\n")
mu			<- rowMeans(shuf_rat,na.rm=TRUE)
sigma		<- apply(shuf_rat,1,sd,na.rm=TRUE)
orig_z		<- (orig_rat - mu)/sigma

orig_pct	<-  sapply(1:length(orig_rat), function(i) {
	tmp		<- na.omit(shuf_rat[i,])
	return(sum(tmp>= orig_rat[i])/length(tmp))
})

maxShufs	<- sapply(1:nrow(shuf_rat),function(i) {
	max(which(!is.na(shuf_rat[i,])))
})

qval	<- p.adjust(orig_pct, method="BH")

out <- data.frame(NETWORK=basename(fList), 
	orig_pp=orig[,1], orig_rest=orig[,2],
	ENR=orig_rat, TOTAL_INT=log10(orig[,1]+orig[,2]),
	numPerm=maxShufs,
	shuf_mu=mu, shuf_sigma=sigma,
	Z=orig_z,pctl=orig_pct,Q=qval)

write.table(out,
			file=sprintf("%s/%s.stats.txt", outDir, outPref),
			sep="\t",col=TRUE,row=F,quote=F)

return(out)

}, error=function(ex) {
	print(ex)
}, finally={
		cat("Closing log\n")
		sink(NULL)
})
}
