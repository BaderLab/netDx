# functions to download and update GeneMANIA jar used for network integration

#' wrapper function for getting BiocFileCache associated with netDx package
#' 
#' @return BiocFileCache object associated with netDx
#' @import BiocFileCache
#' @importFrom rappdirs user_cache_dir
.get_cache <- function() {
    cache <- rappdirs::user_cache_dir(appname = "netDx")
    BiocFileCache::BiocFileCache(cache,ask=FALSE)
}

#' download and update GeneMANIA jar file 
#' 
#' @param verbose (logical) print messages
#' @examples getGMjar_path()
#' @return (char) Path to local cached copy of GeneMANIA jar file..
#' or initial download is required 
#' @export
getGMjar_path <- function(verbose = FALSE) {

	java_ver <- suppressWarnings(
		system2("java", args="--version",stdout=TRUE,stderr=NULL)
	)
	if (any(grep(" 11",java_ver)) || any(grep(" 12",java_ver)) || any(grep(" 13",java_ver)) || any(grep(" 14",java_ver)) || any(grep(" 16",java_ver))) {
		if (verbose) message("Java 11+ detected")
    	fileURL <- paste("https://download.baderlab.org/netDx/java11/", 
			"genemania-netdx.jar",sep="")
	} else {
		if (verbose) message("Java 8 detected")
    	fileURL <- paste("https://download.baderlab.org/netDx/java8/", 
			"genemania-netdx.jar",sep="")
	}
	
    bfc <- .get_cache()
    bfcrpath(bfc, fileURL)
}

#' fetch pathway definitions from downloads.baderlab.org
#' 
#' @details Fetches genesets compiled from multiple curated pathway
#' databases. Downloaded from: http://download.baderlab.org/EM_Genesets/
#' The file contains pathways from HumanCyc, NetPath, Reactome, NCI
#' Curated Pathways and mSigDB.
#' For details see Merico D, Isserlin R, Stueker O, Emili A and GD Bader.
#' (2010). PLoS One. 5(11):e13984.
#' @param verbose (logical) print messages
#' @examples fetchPathwayDefinitions("October",2020)
#' @param day (integer)
#' @param month (numeric or char) month of pathway definition file. Can be
#' numeric or text (e.g. "January","April"). If NULL, fails.
#' @param year (numeric) year of pathway definition file. Must be in
#' yyyy format (e.g. 2018). If NULL, fails.
#' @return (char) Path to local cached copy of GMT file
#' or initial download is required 
#' @importFrom httr HEAD
#' @export
#' @examples 
#' fetchPathwayDefinitions("January",2018)
#' fetchPathwayDefinitions(month=10,year=2020)
fetchPathwayDefinitions <- function(month=NULL,year=NULL,day=1,verbose=FALSE){
	if (is.null(month) || is.null(year)) {
		stop("Please provide a month and year.")
		#month <- month.name[as.integer(format(Sys.Date(),"%m"))]
		#year <- as.integer(format(Sys.Date(),"%Y"))
	}
	if (class(month) %in% c("numeric","integer")) {
		month <- month.name[month]
	}
		pdate <- sprintf("%s_%02d_%i",month,day,year)
    	pathwayURL <- paste("http://download.baderlab.org/EM_Genesets/", 
		sprintf("%s/Human/symbol/",pdate),
        sprintf("Human_AllPathways_%s_symbol.gmt",pdate),
		 sep = "")

	message(sprintf("Fetching %s",pathwayURL))
   	bfc <- .get_cache()
	chk <- httr::HEAD(pathwayURL)
	if (chk$status_code==404) {
		stop(paste(sprintf("The pathway file for %02d %s %i doesn't exist.",day,month,year),
				"Select a different date. ",
				"See http://download.baderlab.org/EM_Genesets/Human/symbol for options.",
				sep=" "))
	}
    bfcrpath(bfc, pathwayURL)
}
