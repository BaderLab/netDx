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
	cat(sprintf("%s\n",orgName), file=outFile)			# org name
	cat(sprintf("%s\n",
				paste(qSamps,collapse="\t")),			# sample names
				file=outFile,append=TRUE)
	cat(sprintf("%s\n",paste(incNets,collapse="\t")),	# networks
				file=outFile,append=TRUE) 				# num2return
	cat(sprintf("%i\n",numReturn),file=outFile,append=TRUE)
	cat("automatic\n",file=outFile,append=TRUE)			# combining
}
