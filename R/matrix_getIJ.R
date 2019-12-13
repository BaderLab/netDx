#' Converts matrix index (1 to m*n) to row (m) and column (n) number
#'
#' @param dimMat (integer vector of length 2) output of \code{dim()} for
#' matrix in question
#' @param idx (integer vector of length n) matrix indices
#' @return (matrix) n-by-2, first column has row indices ; second column
#' has col indices
matrix_getIJ <- function(dimMat, idx) {
    nr <- dimMat[1]
    nc <- dimMat[2]
    
    out <- matrix(NA, nrow = length(idx), ncol = 2)
    out[, 1] <- idx%%nr
    if (any(out[, 1] %in% 0)) {
        out[which(out[, 1] %in% 0)] <- nr
    }
    
    out[, 2] <- ceiling(idx/nr)
    
    out
}
