#' Get ENR for all networks in a specified directory 
#' 
#' @details For each network, compute the number of (+,+) and other 
#' {(+,-),(-,+),(-,-)} interactions. 
#'  From this compute network ENR.
#' The measure of (+,+)-enrichment is defined as: 
#' ENR(network N) = ((num (+,+) edges) - (num other edges))/(num edges).
#' A network with only (+,+) interactions has an ENR=1 ; a network with
#' no (+,+) interactions has an ENR=-1; a network with a balance of the two
#' has ENR=0.
#' @param netDir (char) directory containing interaction networks
#' @param pheno_DF (data.frame) table with patient ID and status.
#'	Must contain columns for Patient ID (named "ID") and class
#' (named "STATUS"). Status should be a char; value of predictor class 
#' should be specified in \code{predClass} param; 
#'	all other values are considered non-predictor class
#' Rows with duplicate IDs will be excluded.
#' @param predClass (char) value for patients in predictor class
#' @param netGrep (char) pattern for grep-ing network text files, used in
#' dir(pattern=..) argument
#' @param enrType (char) how enrichment should be computed. Options are:
#' 1) binary: Skew of number of (+,+) interactions relative to other
#' interactions. Used when all edges in network are set to 1 (e.g. 
#' shared CNV overlap)
#' 2) corr: 0.5*((mean weight of (+,+) edges)-(mean weight of other edges))
#' @param ... arguments for \code{countIntType_batch}
#' @return (list):
#' 1) plusID (char) vector of + nodes
#' 2) minusID (char) vector of - nodes
#' 3) orig_rat (numeric) \code{ENR} for data networks
#' 4) fList (char) set of networks processed
#' 5) orig (data.frame) output of \code{countIntType_batch} for input
#' networks
#' @export
#' @examples
#' data(npheno)
#' netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
#' x <- getEnr(netDir, npheno, "case",netGrep=".txt$")
#' print(x$orig_rat) 
#'
getEnr	<- function(netDir, pheno_DF,predClass,netGrep="_cont.txt$",
	enrType="binary",...) {
if (missing(predClass)) stop("predClass must be supplied.\n")

fList	<- dir(path=netDir,pattern=netGrep)
fList	<- paste(netDir,fList,sep="/")
cat(sprintf("Got %i networks\n", length(fList)))

# get + and - IDs
pheno	<- pheno_DF
pheno	<- pheno[!duplicated(pheno$ID),]
plus_idx	<- which(pheno$STATUS %in% predClass)
plusID	<- pheno$ID[plus_idx]
minusID	<- pheno$ID[setdiff(1:nrow(pheno),plus_idx)]

cat(sprintf("Total %i subjects ; %i of class %s, %i other\n", 
			nrow(pheno), length(plusID), predClass, length(minusID)))

cat("* Computing real (+,+) (+,-)\n")
	# first get original 
	t0 <- system.time(orig <- countIntType_batch(fList,plusID,minusID,
							enrType=enrType,...))
	#print(head(orig))
	print(t0)
	cat("\n")
	if (enrType=="binary") {
		orig_rat	<- (orig[,1]-orig[,2])/(orig[,1]+orig[,2])
	} else if (enrType=="corr") {
		# divide by 2 to rescale to -1 to 1
		orig_rat <- (orig[,1]-orig[,2])/2
	} else { # invalid type
		orig_rat <- NA	
	}
	names(orig_rat)	<- basename(fList)
	
	return(list(plusID=plusID, minusID=minusID,orig_rat=orig_rat,
				fList=fList,orig=orig))
}
