#' platform-specific file separator
#' 
#' @description Returns OS-specific file separator
#' @return (char) "\\" if Windows, else "/"
#' @examples
#' getFileSep()
#' @export
getFileSep <- function(){
  if (.Platform$OS.type=="windows") return("\\")
  else return(.Platform$file.sep)
}

#' Generate random alphanumerical string of length 10
#'
#' @details Used to create multiple temporary directories during an R session
#' @param numStrings (integer) number of strings to generate
#' @return vector of length n, each with 10-char alphanumerical strings
#' @examples
#' randAlphanumString()
#' @export
randAlphanumString <- function(numStrings = 1L) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, numStrings, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, numStrings, TRUE)), 
		sample(LETTERS, numStrings, TRUE))
}
