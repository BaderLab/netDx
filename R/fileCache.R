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
	if (any(grep(" 11",java_ver)) || any(grep(" 12",java_ver)) || any(grep(" 13",java_ver)) || any(grep(" 14",java_ver))) {
		if (verbose) message("Java 11+ detected")
    	fileURL <- paste("http://download.baderlab.org/netDx/java11/", 
			"genemania-netdx.jar",sep="")
	} else {
		if (verbose) message("Java 8 detected")
    	fileURL <- paste("http://download.baderlab.org/netDx/java8/", 
			"genemania-netdx.jar",sep="")
	}
	
    bfc <- .get_cache()
    rid_rec <- bfcquery(bfc, "GM_jar", "rname")
    rid <- rid_rec$rid
	urlChanged <- FALSE
    if (!length(rid)) {
        if (verbose) 
            message("Downloading GeneMANIA jar file (only required once)")
        rid <- names(bfcadd(bfc, "GM_jar", fileURL))
    } else {
		urlChanged <- rid_rec$fpath != fileURL
	}
	if (urlChanged) {
		message("File path changed. Updating...")
		bfcupdate(bfc, rid, fpath=fileURL, ask=FALSE)
	}
    if (!isFALSE(bfcneedsupdate(bfc, rid))){
		bfcdownload(bfc, rid,ask=FALSE)
	}
    
    bfcrpath(bfc, rids = rid)
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
#' @examples fetchPathwayDefinitions()
#' @param month (char) month of pathway definition file. Must be
#' textual name (e.g. "January","April"). If NULL, gets the latest
#' release
#' @param year (numeric) year of pathway definition file. Must be in
#' yyyy format (e.g. 2018).
#' @return (char) Path to local cached copy of GMT file
#' or initial download is required 
#' @export
#' @examples 
#' fetchPathwayDefinitions("January",2018)
#' fetchPathwayDefinitions()
fetchPathwayDefinitions <- function(month=NULL,year=NULL,verbose=FALSE){
	if (is.null(month) || is.null(year)) {
		month <- month.name[as.integer(format(Sys.Date(),"%m"))]
		year <- as.integer(format(Sys.Date(),"%Y"))
	}
		pdate <- sprintf("%s_01_%i",month,year)
    pathwayURL <- paste("http://download.baderlab.org/EM_Genesets/", 
		sprintf("%s/Human/symbol/",pdate),
        sprintf("Human_AllPathways_%s_symbol.gmt",pdate),
		 sep = "")

	message(sprintf("Fetching %s",pathwayURL))
    bfc <- .get_cache()
	cache_name <- sprintf("%s_Pathways",pdate)
    rid <- bfcquery(bfc, cache_name,  "rname")$rid
    if (!length(rid)) {
        if (verbose) 
            message(paste("Downloading example pathway definitions (only ", 
								"required once)", 
                sep = ""))
        rid <- names(bfcadd(bfc, cache_name, pathwayURL))
    }

    if (!isFALSE(bfcneedsupdate(bfc, rid))) 
        bfcdownload(bfc, rid)
    
    bfcrpath(bfc, rids = rid)
}
