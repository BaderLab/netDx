#' Similarity metric of normalized difference
#'
#' @details Similarity metric used when data for a network consists of
#' exactly 1 continuous variable  (e.g. a network based only on "age").
#' When number of variables is 2-5, use avgNormDiff() which
#' takes the average of normalized difference for individual variables
#' @param (numeric) vector of values, one per patient (e.g. ages)
normDiff <- function(x) {
    #if (nrow(x)>=1) x <- x[1,]
    nm <- colnames(x)
    x <- as.numeric(x)
    n <- length(x)
    rngX  <- max(x,na.rm=T)-min(x,na.rm=T)

    out <- matrix(NA,nrow=n,ncol=n);
    # weight between i and j is
    # wt(i,j) = 1 - (abs(x[i]-x[j])/(max(x)-min(x)))
    for (j in 1:n) out[,j] <- 1-(abs((x-x[j])/rngX))
    rownames(out) <- nm; colnames(out)<- nm
    out
}
