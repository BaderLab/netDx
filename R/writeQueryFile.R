#' Wrapper to write GeneMANIA query file
#'
#' @param qSamps (char) vector of patient IDs in query
#' @param incNets (char) vector of networks to include in this analysis
#' (features/pathway names). Useful for subset-based feature selection
#' @param numReturn (integer) number of patients to return in ranking file
#' @param outFile (char) path to output file
#' @param orgName (char) organism name
#' @return No value. Side effect of writing the query file to
#' \code{outFile}
#' @examples
#' data(xpr,pheno,cnv_GR)
#' writeQueryFile(pheno$ID[1:5], "all",nrow(pheno), "myquery.txt")
#' @export
writeQueryFile <- function(qSamps, incNets="all", numReturn=1L, outFile,
  orgName="predictor") {
	fileConn <- file(outFile,"w")
	writeLines(sprintf("%s",orgName),con=fileConn)# org name
	writeLines(sprintf("%s",paste(qSamps,collapse="\t")),con=fileConn)
	writeLines(sprintf("%s",paste(incNets,collapse="\t")),con=fileConn)# networks 	
	writeLines(sprintf("%i",numReturn),con=fileConn) #num2return
	writeLines("automatic",con=fileConn)			# combining
	close(fileConn)
}
