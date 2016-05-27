#' Clean pathway name so it can be a filename.
#' 
#' @param curP (char) pathway name
#' @export
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

	return(pforfile)
}
