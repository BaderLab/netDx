#' evaluate test performance for CNV-based predictor
#'
#' @param testPheno (data.frame) sample metadata table with ID and STATUS
#' @param test_p_GR (GRanges) genomic events in patients. Must have ID column
#' that maps an event to a patient, and LOCUS_NAMES with comma-separated names
#' of units (e.g. genes) that the ranges overlap
#' @param selFeature_GR (list of GRanges) genomic ranges for selected features
#' NOTE: These are the features over which patients will be evaluated so 
#' limit to nets that pass feature selection
#' @param outDir (char) root directory in which data must be saved
#' @export
resampling_predTest_CNV <- function(testPheno, test_p_GR, selFeature_GR,outDir) {
	
netDir <- sprintf("%s/networks_test",outDir)

# get overlap of patients with genes in pathway

browser()

}
