#' various similarity functions

#' @param (data.frame) measures by patients
#' @param K, W (doc tba)
#' @export
sim.pearscale <- function (dat, K = 20, alpha = 0.5) {
ztrans <- function(m) {
	m <- as.matrix(m)
	m2 <- apply(m,1,function(x) { (x-mean(x,na.rm=T))/sd(x,na.rm=T)})
	m2
}
normalize <- function(X) {
		print(dim(X))
        row.sum.mdiag <- rowSums(X,na.rm=T) - diag(X)
        row.sum.mdiag[row.sum.mdiag == 0] <- 1
        X <- X/(2 * (row.sum.mdiag))
        diag(X) <- 0.5
		X <- (X+t(X))/2
        return(X)
}

	if (nrow(dat)<6) {	
		z1 <- ztrans(dat)
		euc <- as.matrix(dist(z1,method="euclidean"))^(1/2)
	} else {
		euc <- as.matrix(1-cor(dat,method="pearson",use="pairwise.complete.obs"))
	}
    N <- nrow(euc)
    euc <- (euc + t(euc))/2
    sortedColumns <- as.matrix(t(apply(euc, 2, sort,na.last=TRUE)))
	print(dim(sortedColumns))
    finiteMean <- function(x) {
        return(mean(x[is.finite(x)],na.rm=T))
    }
    means <- apply(sortedColumns[, 1:K + 1], 1, finiteMean) +
        .Machine$double.eps
    avg <- function(x, y) {
        return((x + y)/2)
    }
    Sig <- outer(means, means, avg)/3 * 2 + euc/3 + .Machine$double.eps
    Sig[Sig <= .Machine$double.eps] <- .Machine$double.eps
    densities <- dnorm(euc, 0, alpha * Sig, log = FALSE)

    W <- (densities + t(densities))/2
	W <- normalize(W)
	# remove patients with no datapoints (full column/row of NAs)
	idx <- which(rowSums(is.na(euc))==ncol(W)-1)
	if (any(idx)) { 
		W <- W[-idx,]
		idx <- which(colSums(is.na(euc))==ncol(W)-1)
		W <- W[,-idx]
	}
    return(W)
}
