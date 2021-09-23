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
  return(c("pearson_corr","normDiff","avgNormDiff",
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


makeNetFunc <- function(dataList, groupList, netDir, sims,...){    
    settings <- list(dataList=dataList,groupList=groupList,
                    netDir=netDir,sims=sims)
    netList <- c()    
    for (nm in names(sims)){
        csim <- sims[[nm]]
        netList_cur <- NULL

        cur_set <- settings; 
        cur_set[["name"]] <- nm; cur_set[["similarity"]] <- csim

        if (!is.null(groupList[[nm]])){
            if (class(csim)=="function") {# custom function
                netList_cur <- psn__custom(cur_set,csim,...)
            } else if (csim == "pearson_corr") {
                netList_cur <- psn__corr(cur_set,...)
            } else {
                netList_cur <- psn__builtIn(cur_set,...)
            }
            netList <- c(netList,netList_cur)
        }
    }
    unlist(netList)
}

#' make PSN for built-in similarity functions
#'
#' @param settings (list) from makeNetFunc
psn__builtIn <- function(settings,...){
funcs <- list(
    "normDiff"=normDiff,
    "avgNormDiff"=avgNormDiff,
    "sim.pearscale"=sim.pearscale,
    "sim.eucscale"=sim.eucscale
)

    message(sprintf("Layer %s: Function %s",settings$name,settings$similarity))

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
psn__custom <- function(settings,fn, ...){
    nm <- settings$name
    message(sprintf("Layer %s: CUSTOM FUNCTION",settings$name))
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
psn__corr <- function(settings,...){
    message(sprintf("Layer %s: PEARSON CORR",settings$name))
    nm <- settings$name
    netList <- makePSN_NamedMatrix(
				settings$dataList,
				rownames(settings$dataList[[nm]]),	## names of measures (e.g. genes, CpGs)
				settings$groupList[[nm]],			## how to group measures in that layer
				settings$netDir,						## leave this as-is, netDx will figure out where this is.
				verbose=FALSE, 			
				writeProfiles=TRUE,   		## use Pearson correlation-based similarity
				...
				)
    return(netList)
}