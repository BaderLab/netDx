#' simple capitalization
#' @details used to format feature names so they are not in all-caps
#' @param x (char) name
#' @return (char) Changes case so start of each word is in upper-case, and
#' the rest is in lowercase
#' @examples simpleCap("this IS a TEST sEnTenCe")
#' @export
simpleCap <- function(x) {
	x <- tolower(x)
	s <- strsplit(x, " ")[[1]]
	x <- paste(toupper(substring(s, 1, 1)), substring(s, 2),
	sep = "", collapse = " ")
	x 
}
