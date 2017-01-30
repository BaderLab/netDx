#' Split samples into training and testing set
#'
#' @param pheno_DF (data.frame) patient information
#' Must contain the following columns:
#' 1. ID: (char) patient IDs
#' 2. STATUS: (char) patient classes. Values not equal to \code{predClass}
#' will be considered as "other"
#' Expects rows with unique IDs
#' @param pctT (numeric between 0 and 1) Fraction of patients to randomly
#' assign to the training set. The remainder will be used for blind test 
#' set
#' @param setSeed (integer) Random number generator seed for splitting
#' the training and test set. If not NULL, seed will be set to this value.
#' @param verbose (logical) print messages
#' @return (char) vector of length \code{nrow(pheno_DF)}, with values of 
#' "TRAIN" or "TEST". The order corresponds to pheno_DF; a patient labelled
#' "TRAIN" has been assigned to the training set, and one labelled "TEST"
#' as been assigned to the test set.
#' @examples
#' data(TCGA_mini)
#' x <- splitTestTrain(pheno)
#' @export
splitTestTrain <- function(pheno_DF,pctT=0.7,setSeed=42,verbose=FALSE) { 

if (!is.null(setSeed)) {
	cat(sprintf("Setting seed for reproducibility: %i\n",setSeed))
	set.seed(setSeed); # make reproducible
}

lvls <- unique(pheno_DF$STATUS)
IS_TRAIN	<- rep("TEST",nrow(pheno_DF))
for (lv in lvls) {
	idx <- which(pheno_DF$STATUS %in% lv)
	IS_TRAIN[sample(idx, floor(pctT*length(idx)), F)] <- "TRAIN"
}

###plus_idx	<- which(pheno_DF$STATUS %in% predClass)
###other_idx	<- setdiff(1:nrow(pheno_DF),plus_idx)
#### randomly assign test/train
###IS_TRAIN[sample(other_idx, floor(pctT*length(other_idx)), F)] <- "TRAIN"

IS_TRAIN <- factor(IS_TRAIN,levels=c("TRAIN","TEST"))

pheno_DF	<- cbind(pheno_DF,IS_TRAIN=IS_TRAIN)
print(table(pheno_DF[,c("STATUS","IS_TRAIN")]))


return(IS_TRAIN)
}
