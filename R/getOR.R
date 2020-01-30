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
#' @param netList (char) vector of networks of interest (e.g. those 
#' passing feature selection)
#' @param verbose (logical) print messages
#' @return List. 1) stats: statistics on group overlap with netList,
#' This is a 2xK matrix, where rows are classes (predClass,other), and 
#' columns are: total samples, samples overlapping nets, % overlap
#' 2) relEnr: relative enrichment of \code{predClass} over other
#' @examples
#' data(npheno)
#' netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
#' x <- countPatientsInNet(netDir,dir(netDir,pattern="txt$"), npheno[,1])
#' getOR(x,npheno,"case",colnames(x)[1]) # should give large RelEnr
#' getOR(x,npheno,"case",colnames(x)[2]) # should give RelEnr of 0
#' getOR(x,npheno,"case",colnames(x)[3]) # should give RelEnr of 1
#' @export
getOR <- function(pNetworks, pheno_DF, predClass, netList,verbose=TRUE) {
# TODO this function only makes sense in the context of structural
# variants. Not e.g. in the case of continuous valued data like gene 
# expression. The name should perhaps reflect that.

predSamps	<- pheno_DF$ID[pheno_DF$STATUS %in% predClass]
otherSamps	<- pheno_DF$ID[!pheno_DF$STATUS %in% predClass]

#print(length(netList))
# limit universe to netList
idx			<- which(colnames(pNetworks)%in% netList)
if (length(idx)<1) 
		return(out=list(stats=matrix(NA,nrow=2,ncol=3),relEnr=NA,
						OLsamps=NA))

pNetworks 	<- pNetworks[,idx,drop=FALSE]
# limit patients to those that overlap netList
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
	cat(sprintf("Relative enrichment of %s: %1.3f\n", predClass, relEnr))

out <- list(stats=outmat, relEnr=relEnr,OLsamps=OLsamps)

out
}
