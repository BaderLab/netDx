#' Clean pathway name so it can be a filename.
#' 
#' @param curP (char) pathway name
#' @export
#' @examples
#' cleanPathwayName("7-(3-AMINO-3-CARBOXYPROPYL)-WYOSINE BIOSYNTHESIS%HUMANCYC%PWY-7286") 
#' @return (char) Cleaned pathway name
cleanPathwayName <- function(curP) {
	pforfile	<- gsub(" ","_",curP)
	pforfile	<- gsub("<","_",pforfile)
	pforfile	<- gsub(">","_",pforfile)
	pforfile	<- gsub("\\(","_",pforfile)
	pforfile	<- gsub("\\)","_",pforfile)
	pforfile	<- gsub("&","_",pforfile)
	pforfile	<- gsub(";","_",pforfile)
	pforfile	<- gsub("\\/","_",pforfile)
	pforfile	<- gsub("\\\354","X",pforfile)
	pforfile	<- gsub("\\\302\\\240","_",pforfile)
	pforfile	<- gsub("\\\240","X",pforfile)
	pforfile	<- gsub("\\+","plus",pforfile)

	return(pforfile)
}
