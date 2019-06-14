#' evaluate test performance for CNV-based predictor
#'
#' @param testPheno (data.frame) sample metadata table with ID and STATUS
#' @param test_p_GR (GRanges) genomic events in patients. Must have ID column
#' that maps an event to a patient, and LOCUS_NAMES with comma-separated names
#' of units (e.g. genes) that the ranges overlap
#' @param selFeature_GR (GRanges) genomic ranges for selected features
#' NOTE: These are the features over which patients will be evaluated so 
#' limit to nets that pass feature selection. An example is the set of genes
#' belonging to feature-selected pathways
#' @param predClass (char) label for positve class
#' @param cliqueNetGenes (GRanges) optional. If provided, universe is 
#' limited to patients overlapping ranges in this set. e.g. the set of genes
#' in all pathways that pass feature selection
#' @param verbose (logical) print messages
#' @return (list) 
#' 1) pos (char): vector of patients called +
#' 2) neg (char): vector of patients called -
#' 3-6) tp,tn,acc,ppv (numeric): num true pos, num true neg, accuracy and
#' positive predictive value
#' @export
resampling_predTest_CNV <- function(testPheno, test_p_GR, selFeature_GR,
	predClass,cliqueNetGenes=NULL,verbose=TRUE) {

if (!is.null(cliqueNetGenes)) {
	cat("*** Limiting universe ***\n")
	test_p_GR <- subsetByOverlaps(test_p_GR, cliqueNetGenes)
	id <- unique(test_p_GR$ID)
	cat(sprintf("%i of %i patients left\n", length(id), nrow(testPheno)))
	testPheno <- subset(testPheno,ID %in% id)
	
}

# get overlap of patients with genes in pathway
has_ol <- rep(0,nrow(testPheno))
names(has_ol) <- testPheno$ID
ol <- subsetByOverlaps(test_p_GR, selFeature_GR)
has_ol[which(names(has_ol) %in% ol$ID)] <- 1

called_pos <- names(has_ol)[which(has_ol>0)]
called_neg <- names(has_ol)[which(has_ol<1)]
tp <- intersect(called_pos, 
	testPheno$ID[which(testPheno$STATUS %in% predClass)])
tp <- length(tp)

tn <- intersect(called_neg, 
	testPheno$ID[which(!testPheno$STATUS %in% predClass)])
tn <- length(tn)

acc <- (tp+tn)/nrow(testPheno)
ppv <- tp/length(called_pos)
real_pos <- sum(testPheno$STATUS %in% predClass)
real_neg <- sum(!testPheno$STATUS %in% predClass)
tpr <- tp/real_pos
fpr <- (length(called_pos)-tp)/real_neg


out <- list(called_pos=called_pos,called_neg=called_neg,
			pos=real_pos,neg=real_neg,
			tp=tp,tn=tn,tpr=tpr,fpr=fpr,
			acc=acc,ppv=ppv)
out
}
