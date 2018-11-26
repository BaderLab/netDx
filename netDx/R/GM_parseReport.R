#' parses GeneMANIA query report into the PRANK and NRANK tables
#'
#' @param resFile (char) path to GeneMANIA query results
#' @return No value. Side effect of writing PRANK and NRANK files to the 
#' same directory as resFile. Suffixes filename with PRANK and NRANK 
#' respectively
#' @export
#' @examples
#' GM_db <- sprintf("%s/extdata/GM_db", path.package("netDx"))
#' GM_query <- sprintf("%s/extdata/GM_query.txt",path.package("netDx"))
#' x <- runGeneMANIA(GM_db,GM_query,"/tmp")
#' GM_parseReport(x)
#' 
GM_parseReport <- function(resFile) {
	os <- Sys.info()['sysname']
	# use OS/X version of tail
	if (any(grep("Darwin",os,ignore.case=TRUE))) {
	system(sprintf(
		"sed -n '/Score/,/Network/p' %s | tail -r | tail -n+3 | tail -r > %s.PRANK",
		   resFile,resFile))
	system(sprintf(
		"sed -n '/^Network/,/Gene/p' %s | sed '/^Gene/q' | tail -r | tail -n+3 | tail -r > %s.NRANK",
		resFile,resFile));

	# use *nix version of tail
	} else if (any(grep("n[iu]x",os,ignore.case=TRUE))){
    # Looking for column names "Network Group	Network	Weight	Title	Authors	Year	Publication	PMID	URL	Processing Method	Interactions	Source	Source URL	Tags"
    # only until occurence of "Gene 1	Gene 2	Weight	Type	Source" header - 
    # indicating gene interactions
    # later we further crop away "dummy_group		100.00" in networkTally for NRANK
  	system(sprintf("sed -n '/^Network/,/Gene/p' %s | sed '/^Gene/q' | head -n -2 | tail -n+1> %s.NRANK",
  		resFile,resFile));
  	
  	system(sprintf("sed -n '/Gene\tScore/,/Network/p' %s | head -n -2 > %s.PRANK",
  		resFile,resFile))
	} else {
		stop("Extracting NRANK/PRANK doesn't currently work on Win")
	}
}
