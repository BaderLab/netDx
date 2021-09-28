#' Assign train/test labels over several resamplings of the data.
#' 
#' @details This function is useful when feature selection needs to 
#' occur over multiple resamplings of the data, as a strategy to reduce
#' overfitting. Each sample serves as a test for exactly one resampilng, 
#' and as a training sample for the others. The method is provided with the
#' positive label and splits the samples so that an even number of positive
#' and negative classes are represented in all the resamplings (i.e. it
#' avoids the situation where one resampling has too many positives and 
#' another has too few).
#' @param pheno_DF (data.frame) table with patient ID and status.
#' Must contain columns for Patient ID (named 'ID') and class
#' (named 'STATUS'). Status should be a char; value of predictor class 
#' should be specified in \code{predClass} param; 
#' all other values are considered non-predictor class
#' Expects rows with unique IDs
#' Rows with duplicate IDs will be excluded.
#' @param nFold (integer) number of resamplings. Each sample will be a test
#' sample exactly once.
#' @param predClass (char) name of predictor class
#' @param verbose (logical) print messages
#' @return (list) of length nFold, each with char vector of length 
#' nrow(pheno_DF). Values of 'TRAIN' or 'TEST'
#' @examples
#' data(pheno) 
#' x <- splitTestTrain_resampling(pheno,predClass='LumA')
#' @export
splitTestTrain_resampling <- function(pheno_DF, nFold = 3L, predClass, 
		verbose = FALSE) {
    
    plus_idx <- which(pheno_DF$STATUS %in% predClass)
    other_idx <- setdiff(seq_len(nrow(pheno_DF)), plus_idx)
    
    # num +/- that should be test per resampling
    plus_csize <- floor((1/nFold) * length(plus_idx))
    other_csize <- floor((1/nFold) * length(other_idx))
    
    plus_tsize <- length(plus_idx) - plus_csize
    other_tsize <- length(other_idx) - other_csize
    
    if (verbose) {
        message(sprintf("\t(+) %s : %i total ; %i train, %i held-out per\n", 
						predClass, length(plus_idx), plus_tsize, plus_csize))
        message(sprintf("\t(-) (!%s): %i total ; %i train, %i held-out per\n", 
						predClass, length(other_idx), other_tsize, other_csize))
    }
    
    # randomize order for test assignment
    plus_order <- sample(plus_idx, replace = FALSE)
    other_order <- sample(other_idx, replace = FALSE)
    
    out <- list()
    for (k in seq_len(nFold)) {
        status <- rep("TRAIN", nrow(pheno_DF))
        
        # first for + samples
        sidx <- ((k - 1) * plus_csize) + 1
        eidx <- k * plus_csize
        if (k == nFold) 
            eidx <- length(plus_idx)
        
        if (verbose) 
            message(sprintf("\t%i (+): %i test (%i-%i);\n", k, 
								eidx - sidx + 1, sidx, eidx))
        status[plus_order[sidx:eidx]] <- "TEST"
        
        # then for - samples
        sidx <- ((k - 1) * other_csize) + 1
        eidx <- k * other_csize
        if (k == nFold) 
            eidx <- length(other_idx)
        
        if (verbose) 
            message(sprintf("\t\t%i (-): %i test\n", k, eidx - sidx + 1))
        status[other_order[sidx:eidx]] <- "TEST"
        
        out[[k]] <- status
    }
    
    out
}


#' Subsample a hold-out set from a larger patient dataset
#'
#' @details Creates a partition of data to be used for model validation after initial model building.
#' In netDx, buildPredictor() is used for model training, and selected features from this exercise
#' are used to validate a held-out dataset with the predict() function. Note that this function
#' identifies a random subsample, which may result in a validation sample that is not representative
#' of your training bias. If this method is used, please use data exploration techniques (e.g. UMAP)
#' to ensure that validation accuracy is not confounded by stratification.
#' @param dataMAE (MultiAssayExperiment) patient data to be subsampled. Must have columns ID (patient ID)
#' and STATUS
#' @param pctValidation (numeric) Fraction of dataset to include in the validation set. Value from 0.05 to 0.95.
#' @param verbose (logical)
#' @return (list) Keys are trainMAE and validationMAE. These contain corresponding MultiAssayExperiments for
#' training and test data
#' @export
subsampleValidationData <- function(dataMAE,pctValidation=0.2,verbose=TRUE) {
    if (pctValidation < 0.05 || pctValidation > 0.95) stop("pctValidation should be between 0.05 and 0.95.")
    if (missing(dataMAE)) stop("Supply dataMAE.")
    if (class(dataMAE)!= "MultiAssayExperiment") stop("dataMAE must be an object of type MultiAssayExperiment.")

    pheno <- colData(dataMAE)
    st <- unique(pheno$STATUS)
    nsamp <- round(pctValidation/length(st) * nrow(pheno))

    idx_holdout <- c()
    for (k in unique(pheno$STATUS)) {
        idx_holdout <- c(idx_holdout, 
            sample(which(pheno$STATUS == k),nsamp,FALSE)
        )  
    }

    holdout <- dataMAE[,rownames(pheno)[idx_holdout]]
    colData(holdout)$ID <- as.character(colData(holdout)$ID)
    tokeep <- setdiff(1:nrow(pheno),idx_holdout)

    dataMAE <- dataMAE[,rownames(pheno)[tokeep]]

    return(list(trainMAE=dataMAE, validationMAE=holdout))
}