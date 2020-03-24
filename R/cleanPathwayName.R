#' Clean pathway name so it can be a filename.
#' 
#' @param curP (char) pathway name
#' @export
#' @examples
#' cleanPathwayName('7-(3-AMINO-3-CARBOXYPROPYL)-WYOSINE BIOSYNTHESIS%HUMANC')
#' @return (char) Cleaned pathway name
cleanPathwayName <- function(curP) {
    pforfile <- gsub(" ", "_", curP)
    pforfile <- gsub("<", "_", pforfile)
    pforfile <- gsub(">", "_", pforfile)
    pforfile <- gsub("\\(", "_", pforfile)
    pforfile <- gsub("\\)", "_", pforfile)
    pforfile <- gsub("&", "_", pforfile)
    pforfile <- gsub(";", "_", pforfile)
    pforfile <- gsub(":", "_", pforfile)
    pforfile <- gsub("\\/", "_", pforfile)
    pforfile <- gsub("\\\xec", "X", pforfile)
    pforfile <- gsub("\\\xc2\\\xa0", "_", pforfile)
    pforfile <- gsub("\\\xa0", "X", pforfile)
    pforfile <- gsub("\\\xca", "_", pforfile)
    pforfile <- gsub("\\+", "plus", pforfile)
    
    return(pforfile)
}
