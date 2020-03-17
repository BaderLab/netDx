#' Get relative proportion of patient classes that contribute to a set of
#' networks
#'
#' @details Feature selected networks should have the property of being
#' enriched in the class of interest; e.g. be enriched in 'case' relative
#' to 'control'. When given a list of networks N, this method computes the
#' number and proportion of patients that overlap N. A high relative 
#' fraction of the predicted class indicates successful feature selection.
#' To create a ROC or precision-recall curve, several calls can be made
#' to this function, one per cutoff.
#' @param pNetworks (matrix) rows are patients, columns are network file
#' filenames. a[i,j] = 1 if patient i has a structural variant in network
#' j; else a[i,j] = 0
#' @param pheno_DF (data.frame) Column "ID" has unique patient identifiers;
#' column "STATUS" has patient class
#' @param predClass (char) Class for which predictor is being built
#' @param netFile (char) vector of networks of interest (e.g. those 
#' passing feature selection)
#' @param verbose (logical) print messages
#' @return List. 1) stats: statistics on group overlap with ,
#' This is a 2xK matrix, where rows are classes (predClass,other), and 
#' columns are: total samples, samples overlapping nets, % overlap
#' 2) relEnr: relative enrichment of \code{predClass} over other
#' @examples
#' d <- tempdir()
#' options(stringsAsFactors=FALSE)
#' pids <- paste("P",seq_len(5),sep="")
#' pheno <- data.frame(ID=pids,STATUS=c(rep("case",3),rep("control",2)))
#' 
#' # write PSN
#' m1 <- matrix(c("P1","P1","P2","P2","P3","P4",1,1,1),byrow=FALSE,ncol=3)
#' write.table(m1,file=sprintf("%s/net1.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' m2 <- matrix(c("P3","P4",1),nrow=1)
#' write.table(m2,file=sprintf("%s/net2.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#'
#' # compute enrichment
#' x <- countPatientsInNet(d,dir(d,pattern=c("net1.txt","net2.txt")), pids)
#' getOR(x,pheno,"case",colnames(x)) # should give large RelEnr
#' @export
getOR <- function(pNetworks, pheno_DF, predClass, netFile ,verbose=TRUE) {
# TODO this function only makes sense in the context of structural
# variants. Not e.g. in the case of continuous valued data like gene 
# expression. The name should perhaps reflect that.

predSamps	<- pheno_DF$ID[pheno_DF$STATUS %in% predClass]
otherSamps	<- pheno_DF$ID[!pheno_DF$STATUS %in% predClass]

#print(length())
# limit universe to 
idx			<- which(colnames(pNetworks)%in% netFile)
if (length(idx)<1) 
		return(out=list(stats=matrix(NA,nrow=2,ncol=3),relEnr=NA,
						OLsamps=NA))

pNetworks 	<- pNetworks[,idx,drop=FALSE]
# limit patients to those that overlap 
OLsamps		<- rownames(pNetworks)[which(rowSums(pNetworks)>=1)]

OLpred		<- sum(OLsamps %in% predSamps)
OLother		<- sum(OLsamps %in% otherSamps)

pctPred		<- OLpred/length(predSamps)
pctOther	<- OLother/length(otherSamps)

if (pctPred < .Machine$double.eps) pctPred <- .Machine$double.eps
if (pctOther < .Machine$double.eps) pctOther <- .Machine$double.eps

relEnr		<- pctPred/pctOther

outmat <- matrix(nrow=2,ncol=3)
colnames(outmat) <- c("total","num OL","pct OL")
rownames(outmat) <- c(predClass, "(other)")
# cases - total, overlapping selPath, fraction
# ctrl - total, overlapping selPath, fraction
outmat[1,] <- c(length(predSamps), OLpred,round(pctPred*100,digits=1))
outmat[2,] <- c(length(otherSamps),OLother,round(pctOther*100,digits=1))

if (verbose) print(outmat)
if (verbose)
	message(sprintf("Relative enrichment of %s: %1.3f", 
		predClass, relEnr))
out <- list(stats=outmat, relEnr=relEnr,OLsamps=OLsamps)
out
}

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
#' d <- tempdir()
#' options(stringsAsFactors=FALSE)
#' pids <- paste("P",seq_len(5),sep="")
#' pheno <- data.frame(ID=pids,STATUS=c(rep("case",3),rep("control",2)))
#' 
#' # write PSN
#' m1 <- matrix(c("P1","P1","P2","P2","P3","P4",1,1,1),byrow=FALSE,ncol=3)
#' write.table(m1,file=sprintf("%s/net1.nettxt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' m2 <- matrix(c("P3","P4",1),nrow=1)
#' write.table(m2,file=sprintf("%s/net2.nettxt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#'
#' # compute enrichment
#' x <- countPatientsInNet(d,dir(d,pattern=c("net1.nettxt","net2.nettxt")), pids)
#' getEnr(d,pheno,"case","nettxt$")
getEnr	<- function(netDir, pheno_DF,predClass,netGrep="_cont.txt$",
	enrType="binary",...) {
if (missing(predClass)) stop("predClass must be supplied.\n")

fList	<- dir(path=netDir,pattern=netGrep)
fList	<- paste(netDir,fList,sep="/")
message(sprintf("Got %i networks", length(fList)))

# get + and - IDs
pheno	<- pheno_DF
pheno	<- pheno[!duplicated(pheno$ID),]
plus_idx	<- which(pheno$STATUS %in% predClass)
plusID	<- pheno$ID[plus_idx]
minusID	<- pheno$ID[setdiff(seq_len(nrow(pheno)),plus_idx)]

message(sprintf("Total %i subjects ; %i of class %s, %i other", 
			nrow(pheno), length(plusID), predClass, length(minusID)))

message("* Computing real (+,+) (+,-)")
	# first get original 
	t0 <- system.time(orig <- countIntType_batch(fList,plusID,minusID,
							enrType=enrType,...))
	#print(head(orig))
	print(t0)
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

#' Counts the number of (+,+) and (+,-) interactions in a single network
#' 
#' @param inFile (char) path to interaction networks
#' @param plusID (char) vector of + nodes
#' @param minusID (char) vector of - nodes
#' @return (numeric of length 2) Number of (+,+) interactions, and 
#' non-(+,+) interactions
#' (i.e. (+,-) and (-,-) interactions)
#' @export
#' @examples
#' d <- tempdir()
#' # write PSN
#' m1 <- matrix(c("P1","P1","P2","P2","P3","P4",1,1,1),byrow=FALSE,ncol=3)
#' write.table(m1,file=sprintf("%s/net1.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#'
#' countIntType(sprintf("%s/net1.txt",d),c("P1","P2","P3"),
#'	c("P4","P5"))
countIntType <- function(inFile, plusID, minusID) { 
	dat <- read.delim(inFile,sep="\t",header=FALSE,as.is=TRUE)
	pp	<- sum(dat[,1] %in% plusID & dat[,2] %in% plusID)

	return(c(pp,nrow(dat)-pp))
}

#' Counts number of (+,+) and (+,-) interactions in a set of networks
#' 
#' @param inFiles (char) path to interaction networks to process
#' @param plusID (char) IDs of + nodes
#' @param minusID (char) IDs of - nodes
#' @param tmpDir (char) path to dir where temporary files can be stored
#' @param enrType (char) see getEnr.R
#' @param numCores (integer) number of cores for parallel processing
#' @import bigmemory
#' @return (matrix) two columns, one row per network 
#' If \code{enrType="binary"}, number of (+,+) and other interactions
#' Otherwise if \code{enrType="corr"} mean edge weight of (+,+) edges and
#' of other edges
#' @examples
#' d <- tempdir()
#' # write PSN
#' m1 <- matrix(c("P1","P1","P2","P2","P3","P4",1,1,1),byrow=FALSE,ncol=3)
#' write.table(m1,file=sprintf("%s/net1.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' m2 <- matrix(c("P3","P4",1),nrow=1)
#' write.table(m2,file=sprintf("%s/net2.txt",d),sep="\t",
#'	col.names=FALSE,row.names=FALSE,quote=FALSE)
#' 
#' countIntType_batch(paste(d,c("net1.txt","net2.txt"),sep="/"),
#' 	c("P1","P2","P3"),c("P4","P5"))
#' @export
countIntType_batch <- function(inFiles,plusID, minusID,tmpDir=tempdir(),
	   enrType="binary",numCores=1L){
	bkFile <- sprintf("%s/tmp.bk",tmpDir)
	if (file.exists(bkFile)) file.remove(bkFile)
	out <- big.matrix(NA, nrow=length(inFiles),ncol=2,
					  type="double",backingfile="tmp.bk",
					  backingpath=tmpDir,
					  descriptorfile="tmp.desc")
	cl <- makeCluster(numCores,outfile=sprintf("%s/shuffled_log.txt",tmpDir))
	registerDoParallel(cl)

	k <- 0
	foreach (k=seq_len(length(inFiles))) %dopar% {
		m <-attach.big.matrix(
				sprintf("%s/tmp.desc",tmpDir))
		if (enrType == "binary")
			m[k,]	<- countIntType(inFiles[k], plusID,minusID)
		else if (enrType == "corr")
			m[k,]	<- getCorrType(inFiles[k],plusID,minusID)
	}

	stopCluster(cl)
	
	out	<- as.matrix(out)
	unlink(sprintf("%s/tmp.bk",tmpDir))
	unlink(sprintf("%s/tmp.desc",tmpDir))
	return(out)
}
