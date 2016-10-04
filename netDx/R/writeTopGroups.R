#' writes subset of gmt with feature-selected pathways
#'
#' @param groupList (list) names are input nets, values are member entities
#' e.g. pathways and genes
#' @param scores (data.frame) feature names (NAME) and score from predictor
#' building (SCORE)
#' @param cutoff (integer) only features with scores >= cutoff will be 
#' included
#' @param units (char) units that were present in the data. Not all units 
#' in groupList may be in the data. This function ensures that the units
#' written out are limited to those for which data were available
#' e.g. for a pathway, units are the genes that were present in the 
#' dataset. Genes that were not assayed would not be included in the output.
#' To include all genes (not recommended unless you know what you are doing)
#' , set to "*".
#' @param outFile (char) where output should be written.
#' @export
writeTopGroups <- function(groupList, scores, cutoff,units="*",outFile) {
	pnames <- scores[which(scores$SCORE>=cutoff),1]
	cat(sprintf("%i of %i groups selected\n", length(pnames),nrow(scores)))

	system(sprintf("cat /dev/null > %s",outFile))
	for (p in pnames) {
		if (units[1]!="*"){
			g <- intersect(groupList[[p]],units)
			cat(sprintf("%s\t%s\t%s\n",p,p,paste(g,collapse="\t")),
				file=outFile,append=TRUE)
		}
	}

	cat(sprintf("GMT written to %s\n", outFile))
}
