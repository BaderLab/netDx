#' Sparsifies patient similarity network
#' 
#' @details Extremely dense patient networks (num interactions approaching
#' n choose 2, where n is number of patients) arise from full matrix data
#' such as gene expression measures. Such networks may mostly comprise of
#' uninformative low-edge weights. When combined with tens to thousands
#' of other such networks, they will generate out of memory errors upon
#' running a GeneMANIA query for cross validation. For example, a dataset
#' with 500 patients, which uses gene expression data organized by 1000 
#' pathways, will generate 1,000 networks each with 124,750 interactions.
#' 
#' This function provides a simple rule to sparsify a network, by 
#' retaining the top nearest neighbours for each patient in the network. 
#' The rules for sparsification are based on GeneMANIA's rules for the same
#' algorithm and are as follows:
#' for each patient p:
#'  1. keep the top \code{k} edge weights.
#'  2. sort edges by decreasing weight.
#' 	3. if \code{keepTies=FALSE}, keep only the first occurring instance of
#' 	each tied edge.
#'	4. if \code{keepTies=TRUE}, keep ties up to a specified maximum 
#'		fraction of the total number of patients 
#'		\code{(MAX_PCT*numPatients)}, or 600 interactions, 
#' whichever is less.
#' done
#' 
#' @param net (char or data.frame) If of type char, should path to 
#' tab-separated SIF file containing the network; there should be no header
#' row. If of type data.frame, should have three columns: 
#' 1) patient 1 (char), 2) patient 2 (char), 3) edge weight (numeric)
#' @param outFile (char) path to write sparsified network
#' @param k (integer) for k nearest-neighbours
#' @param MAX_INT (integer) upper bound on number of interactions per 
#' patient
#' @param MAX_PCT (numeric 0-1) upper bound on number of interactions 
#' as fraction of the total number of patients. See Details.
#' @param numPatients (integer) number of patients in the network. See 
#' Details.
#' @param keepTies (logical) keep edge ties. See Details 
#' @param verbose (logical) print messages
#' @param cutoff (value between 0 and 1) min edge value to keep
#' @return No value. Writes sparsified matrix to \code{outFile}
#' @export
#' @examples
#' require(reshape2) # for melt()
#' data(xpr,pheno,cnv_GR)
#' x <- reshape2::melt(cor(xpr)) # patient 1, patient 2, edge weight
#' sparsifyNet(x,outFile="tmp.txt",numPatients=40)
sparsifyNet <- function(net,outFile,k=50L,MAX_INT=600L,MAX_PCT=0.02,
		numPatients,keepTies=TRUE,verbose=TRUE,cutoff=0.3){
if (is(net,"data.frame")){
	dat <- net
} else if (is(net,"character")) { 
	netFile <- net
	dat <- read.delim(netFile,sep="\t",as.is=TRUE,header=FALSE)
}

	dat <- dat[order(dat[,1]),]
	idx <- which(dat[,3]<cutoff)

	if (any(idx)) dat[idx,3] <- NA
	dat <- na.omit(dat)
	curPat <- dat[1,1] # initialize
	sidx <- 1; eidx <- NA;
	ctr <- 1
	
	newCt <- 0
	t0 <- Sys.time()
	system2(sprintf("cat /dev/null > %s",outFile))
	while ((ctr < nrow(dat))) {
		nextPat <- dat[ctr+1,1]
	
		if (nextPat != curPat) {
			eidx <- ctr
			if (verbose) message(sprintf("%s: %i-%i:", curPat,sidx,eidx))
			# process cur pat's interactions. write to file
			totalInter <- dat[sidx:eidx,3]
			names(totalInter) <- dat[sidx:eidx,2]
	
			if (!keepTies) {
				totalInter <- totalInter[!duplicated(totalInter)]		
			}
			totalInter <- sort(totalInter,decreasing=TRUE)
			n <- length(totalInter)
	
			tokeep <- max(k,min(round(MAX_PCT*numPatients),600));
			tokeep <- min(k,n);
	
			if (verbose) {
				n1 <- length(totalInter)
				message(sprintf("%i -> %i ",n1, tokeep))
				if (tokeep < n1) message("*** trimmed")
				message("\n")
			}
			
			outInter <- totalInter[1:tokeep]
			df <- data.frame(P1=curPat,P2=names(outInter),x=outInter)
			write.table(df,file=outFile,sep="\t",
						append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
			newCt <- newCt + nrow(df)
	
			curPat <- nextPat;
			sidx <- ctr+1
		}
		ctr <- ctr+1 
	}
	# last patient
	eidx <- nrow(dat)
	# process cur pat's interactions. write to file
	totalInter <- dat[sidx:eidx,3]
	names(totalInter) <- dat[sidx:eidx,2]
	n1 <- length(totalInter)
	
	if (!keepTies) {
	totalInter <- totalInter[!duplicated(totalInter)]		
	}
	totalInter <- sort(totalInter,decreasing=TRUE)
	n <- length(totalInter)
	
	tokeep <- max(k,min(MAX_PCT*numPatients,600));
	tokeep <- min(k,n);
	
	if (verbose) message(sprintf("%i -> %i\n",n1, tokeep))
	
	outInter <- totalInter[1:tokeep]
df <- data.frame(P1=curPat,P2=names(outInter),x=outInter)
write.table(df,file=outFile,sep="\t",append=TRUE,col=FALSE,row=FALSE,quote=FALSE)
oldCt <- nrow(dat)

message(sprintf("Interactions trimmed from %i to %i  (sparse factor= %1.2f%%)\n", 
			oldCt, newCt,(newCt/oldCt)*100))
message("Time taken:\n")
print(Sys.time()-t0)

}
