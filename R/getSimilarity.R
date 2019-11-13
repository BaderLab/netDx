#' Measures of patient similarity
#' 
#' @param x (matrix) matrix for which pairwise patient similarity is to be
#' computed. Expects one column per patient, and one measurement per row.
#' @param type (character) name of similarity measure. Currently supports 
#' Pearson correlation ("pearson") or a custom measure ("custom")
#' @param customFunc (function) custom similarity function. Only used when 
#' \code{type="custom"}. The function takes \code{x} as first argument and 
#' can take additional argument. It should return a symmetric matrix of 
#' pairwise patient similarities.
#' @param ... parameter for customFunc
#' @return symmetric matrix of size N, where N is number of samples
#' @examples
#' data(xpr) 
#' x <- getSimilarity(xpr) # similarity by Pearson corr
#' mySim <- function(x) cor(x,method="kendall")
#' x <- getSimilarity(xpr,customFunc=mySim) # custom similarity
#' @importFrom stats cor
#' @export
getSimilarity <- function(x, type="pearson",customFunc,...) {
	switch(type, 
		pearson=round(cor(na.omit(x),method="pearson"),digits=3),
		custom=customFunc(x,...)
    )
}
