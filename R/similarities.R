#' various similarity functions

#' Similarity function: Pearson correlation followed by exponential scaling
#' 
#' @description Computes Pearson correlation between patients. A scaled 
#' exponential similarity kernel is used to determine edge weight. The 
#' exponential scaling considers the K nearest neighbours, so that 
#' similarities between non-neighbours is set to zero. Alpha is a 
#' hyperparameter that determines decay rate of the exponential. For details
#' see Wang et al. (2014). Nature Methods 11:333. 
#' @param dat (data.frame) Patient data; rows are measures, columns are 
#' patients.
#' @param K (integer) Number of nearest neighbours to consider (K of KNN)
#' @param alpha (numeric) Scaling factor for exponential similarity kernel. 
#' Recommended range between 0.3 and 0.8.
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' @examples
#' data(xpr)
#' sim <- sim.pearscale(xpr)
#' @importFrom stats dist cor sd dnorm
#' @export
sim.pearscale <- function(dat, K = 20, alpha = 0.5) {
    ztrans <- function(m) {
        m <- as.matrix(m)
        m2 <- apply(m, 1, function(x) {
            (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        })
        m2
    }
    normalize <- function(X) {
        print(dim(X))
        row.sum.mdiag <- rowSums(X, na.rm = TRUE) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
        X <- (X + t(X))/2
        return(X)
    }
    
    if (nrow(dat) < 6) {
        z1 <- ztrans(dat)
        euc <- as.matrix(dist(z1, method = "euclidean"))^(1/2)
    } else {
        euc <- as.matrix(1 - cor(dat, method = "pearson", 
					use = "pairwise.complete.obs"))
    }
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort, na.last = TRUE)))
    print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)], na.rm = TRUE))
    }
    means <- apply(sortedColumns[, seq_len(K) + 1], 1, finiteMean) 
		means <- means + .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)
    
    W <- (densities + t(densities))/2
    W <- normalize(W)
    # remove patients with no datapoints (full column/row of NAs)
    idx <- which(rowSums(is.na(euc)) == ncol(W) - 1)
    if (any(idx)) {
        W <- W[-idx, ]
        idx <- which(colSums(is.na(euc)) == ncol(W) - 1)
        W <- W[, -idx]
    }
    return(W)
}

#' Similarity method. Euclidean distance followed by exponential scaling
#'
#' @description Computes Euclidean distance between patients. A scaled 
#' exponential similarity kernel is used to determine edge weight. The 
#' exponential scaling considers the K nearest neighbours, so that 
#' similarities between non-neighbours is set to zero. Alpha is a 
#' hyperparameterthat determines decay rate of the exponential. For details,
#' see Wang et al. (2014). Nature Methods 11:333. 
#' @param dat (data.frame) Patient data; rows are measures, columns are 
#' patients.
#' @param K (integer) Number of nearest neighbours to consider (K of KNN)
#' @param alpha (numeric) Scaling factor for exponential similarity kernel. 
#' Recommended range between 0.3 and 0.8.
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' @examples
#' data(xpr)
#' sim <- sim.eucscale(xpr)
#' @importFrom stats dist dnorm sd
#' @export
sim.eucscale <- function(dat, K = 20, alpha = 0.5) {
    ztrans <- function(m) {
        m <- as.matrix(m)
        m2 <- apply(m, 1, function(x) {
            (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
        })
        m2
    }
    normalize <- function(X) {
        print(dim(X))
        row.sum.mdiag <- rowSums(X, na.rm = TRUE) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
        X <- (X + t(X))/2
        return(X)
    }
    nnodata <- which(abs(colSums(dat, na.rm = TRUE)) < .Machine$double.eps)
    # if (length(nodata)>0) dat[nodata] <- median(dat) # impute median
    z1 <- ztrans(dat)
    euc <- as.matrix(dist(z1, method = "euclidean"))^(1/2)
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort, na.last = TRUE)))
    print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)], na.rm = TRUE))
    }
    means <- apply(sortedColumns[, seq_len(K) + 1], 1, finiteMean)
    means <- means + .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)
    
    W <- (densities + t(densities))/2
    W <- normalize(W)
    # remove patients with no datapoints (full column/row of NAs)
    idx <- which(rowSums(is.na(euc)) == ncol(W) - 1)
    if (any(idx)) {
        W <- W[-idx, ]
        idx <- which(colSums(is.na(euc)) == ncol(W) - 1)
        W <- W[, -idx]
    }
    return(W)
}

#' Similarity metric of normalized difference 
#'
#' @details Similarity metric used when data for a network consists of
#' exactly 1 continuous variable  (e.g. a network based only on 'age'). 
#' When number of variables is 2-5, use avgNormDiff() which 
#' takes the average of normalized difference for individual variables
#' @param x (numeric) vector of values, one per patient (e.g. ages)
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' @examples
#' sim <- normDiff(rnorm(10))
#' @export
normDiff <- function(x) {
    # if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
    
    out <- matrix(NA, nrow = n, ncol = n)
    # weight between i and j is wt(i,j) = 1 - (abs(x[i]-x[j])/(max(x)-min(x)))
    for (j in seq_len(n)) out[, j] <- 1 - (abs((x - x[j])/rngX))
    rownames(out) <- nm
    colnames(out) <- nm
    out
}

#' takes average of normdiff of each row in x
#' 
#' @param x (numeric) matrix of values, one column per patient (e.g. ages)
#' @return symmetric matrix of size ncol(dat) (number of patients) containing
#' pairwise patient similarities
#' @examples
#' data(xpr)
#' sim <- avgNormDiff(xpr[,seq_len(2)])
#' @export
avgNormDiff <- function(x) {
    # normalized difference x is vector of values, one per patient (e.g. ages)
    normDiff <- function(x) {
        # if (nrow(x)>=1) x <- x[1,]
        nm <- colnames(x)
        x <- as.numeric(x)
        n <- length(x)
        rngX <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
        
        out <- matrix(NA, nrow = n, ncol = n)
        # weight between i and j is wt(i,j) wt(i,j) = 1 -
        # (abs(x[i]-x[j])/(max(x)-min(x)))
        for (j in seq_len(n)) out[, j] <- 1 - (abs((x - x[j])/rngX))
        rownames(out) <- nm
        colnames(out) <- nm
        out
    }
    
    sim <- matrix(0, nrow = ncol(x), ncol = ncol(x))
    for (k in seq_len(nrow(x))) {
        tmp <- normDiff(x[k, , drop = FALSE])
        sim <- sim + tmp
        rownames(sim) <- rownames(tmp)
        colnames(sim) <- colnames(tmp)
    }
    sim <- sim/nrow(x)
    sim
}

#' built-in similarity functions
#'
allowedSims <- function(){
  return(c("pearsonCorr","normDiff","avgNormDiff",
        "sim.pearscale","sim.eucscale"))
}



#' checks if provided similarity functions are valid. Returns error if not
#'
#' @param sims (list) keys are layer names, values are functions or characters (names of built-in similarity functions)
#' @return TRUE if all pass check. Else throws error.
checkSimValid <- function(sims){
    allowed <- allowedSims()
    for (k in names(sims)){    
        if (class(sims[[k]])!="function"){
            if (class(sims[[k]])!="character"){
                stop(paste("Invalid sims datatype. ",
                    "sims entries must be functions or keywords (characters) ",
                    "for built-in similarity functions.",sep=""))
            } else {
                if (!sims[[k]] %in% allowed){
                    stop(paste(
                            sprintf("sims[[%s]] has invalid similarity type:",k),
                            sims[[k]],". ",
                            "Allowed values are: {%s}",
                            paste(allowed,collapse=",")))
                }
            }
        }
    }
    return(TRUE)
}

#' internal test function to check validity of makeNetFunc and sims
#'
#' @details User must provide either makeNetFunc or sims. This function
#' confirms this.
#' @param makeNetFunc (function) makeNetFunc from buildPredictor()
#' @param sims (list) sims from buildPredictor()
#' @param groupList (list) groupList from buildPredictor()s
#' @return (list) cleaned values for makeNetFunc and Sims
checkMakeNetFuncSims <- function(makeNetFunc,sims,groupList){
    if (is.null(makeNetFunc) && is.null(sims)) {
	stop("Provide either makeNetFunc or sims (preferred).")
} 
if (!is.null(makeNetFunc) && !is.null(sims)){
	stop("Provide either makeNetFunc or sims (preferred).")
}

if (!is.null(sims))	{
	if (class(sims)!="list") stop("sims must be a list.")
	if (all.equal(sort(names(sims)),sort(names(groupList)))!=TRUE) 
		stop("names(sims) must match names(groupList).")
}
return(TRUE)
}

#' Create PSN from provided similarities
#'
#' @details Called by CreatePSN_MultiData(), this is the function that converts user-provided
#' simlarity metrics to internal netDx function calls to generate nets.
#' @param dataList (list) patient data, output of dataList2List()
#' @param groupList (list) measure groupings. Keys match assays(dataList) and are usually different data sources. Values for each are a list of 
#' networks with user-provided groupings. See groupList in buildPredictor() for details.
#' @param netDir (char) path to directory where networks are to be created
#' @param sims (list) keys must be identical to those of groupList. Values are either of type character, used for built-in similarity functions, 
#' or are functions, when a custom function is provided.
#' @param verbose (logical) print messages
#' @param ... values to be passed to PSN creation functions such as makePSN_NamedMatrix().
#' @export
createNetFuncFromSimList <- function(dataList, groupList, netDir, sims,
    verbose=TRUE,...){    
    
    if (length(groupList)!= length(sims)){
        stop("groupList and sims need to be of same length.")
    }
    if (all.equal(sort(names(groupList)),sort(names(sims)))!=TRUE){
        stop("names(groupList) needs to match names(sims).")
    }
    settings <- list(dataList=dataList,groupList=groupList,
                    netDir=netDir,sims=sims)

    if (verbose) message("Making nets from sims")
    netList <- c()    
    for (nm in names(sims)){
        csim <- sims[[nm]]
        netList_cur <- NULL
        if (verbose) message(sprintf("\t%s",nm))

        cur_set <- settings; 
        cur_set[["name"]] <- nm; cur_set[["similarity"]] <- csim

        if (!is.null(groupList[[nm]])){
            if (class(csim)=="function") {# custom function
    
                netList_cur <- psn__custom(cur_set,csim, verbose,...)
            } else if (csim == "pearsonCorr") {
                netList_cur <- psn__corr(cur_set,verbose,...)
            } else {
                netList_cur <- psn__builtIn(cur_set,verbose,...)
            }
            netList <- c(netList,netList_cur)
        }
    }
    if (verbose) {
        message("Net construction complete!")
    }
    unlist(netList)
}

#' make PSN for built-in similarity functions
#'
#' @param settings (list) from makeNetFunc
#' @param verbose (logical) print messages
#' @param ... parameters for makePSN_NamedMatrix()
#' @return (char) names of networks created. Side effect of network creation.
psn__builtIn <- function(settings,verbose,...){

funcs <- list(
    "normDiff"=normDiff,
    "avgNormDiff"=avgNormDiff,
    "sim.pearscale"=sim.pearscale,
    "sim.eucscale"=sim.eucscale
)

    if (verbose) message(sprintf("Layer %s: Built-in function %s",
            settings$name,settings$similarity))

    nm <- settings$name
    netList <- makePSN_NamedMatrix(
        settings$dataList[[nm]],
		rownames(settings$dataList[[nm]]),
		settings$groupList[[nm]],
        settings$netDir,
		simMetric="custom",
        customFunc=funcs[[settings$similarity]], # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,...
    )
    netList
}

#' make PSN for custom similarity functions
#'
#' @param settings (list) from makeNetFunc
#' @param fn (function) custom similarity function
#' @param verbose (logical) print messages
#' @param ... parameters for makePSN_NamedMatrix()
#' @return (char) names of networks created. Side effect of network creation.
psn__custom <- function(settings,fn,verbose, ...){
    nm <- settings$name
    if (verbose) message(sprintf("Layer %s: CUSTOM FUNCTION",settings$name))
    netList <- makePSN_NamedMatrix(
        settings$dataList[[nm]],
		rownames(settings$dataList[[nm]]),
		settings$groupList[[nm]],
        settings$netDir,
		simMetric="custom",customFunc=fn, # custom function
		writeProfiles=FALSE,
		sparsify=TRUE,...
    )
    netList
}

#' wrapper for PSNs using Pearson correlation
#'
#' @param settings (list) from makeNetFunc
#' @param verbose (logical) print messages
#' @param ... parameters for makePSN_NamedMatrix()
#' @return (char) names of networks created. Side effect of network creation.
psn__corr <- function(settings,verbose,...){
    if (verbose) message(sprintf("Layer %s: PEARSON CORR",settings$name))
    nm <- settings$name
    netList <- makePSN_NamedMatrix(
				xpr=settings$dataList[[nm]],
				nm=rownames(settings$dataList[[nm]]),
				namedSets=settings$groupList[[nm]],	
				outDir=settings$netDir,	
				verbose=FALSE, 			
				writeProfiles=TRUE,  
				...
				)
    return(netList)
}