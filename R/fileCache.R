# functions to download and update GeneMANIA jar used for network integration

#' wrapper function for getting BiocFileCache associated with netDx package
#' 
#' @return BiocFileCache object associated with netDx
#' @import BiocFileCache
#' @importFrom rappdirs user_cache_dir
.get_cache <- function() {
    cache <- rappdirs::user_cache_dir(appname = "netDx")
    BiocFileCache::BiocFileCache(cache)
}

#' download and update GeneMANIA jar file 
#' 
#' @param verbose (logical) print messages
#' @examples getGMjar_path()
#' @return (char) Path to local cached copy of GeneMANIA jar file..
#' or initial download is required 
#' @export
getGMjar_path <- function(verbose = FALSE) {
    fileURL <- paste("http://download.baderlab.org/netDx/", "genemania-cytoscape-plugin-3.5.0.jar", 
        sep = "")
    
    bfc <- .get_cache()
    rid <- bfcquery(bfc, "GM_jar", "rname")$rid
    if (!length(rid)) {
        if (verbose) 
            message("Downloading GeneMANIA jar file (only required once)")
        rid <- names(bfcadd(bfc, "GM_jar", fileURL))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid))) 
        bfcdownload(bfc, rid)
    
    bfcrpath(bfc, rids = rid)
}

#' download pathway gmt file used for several examples
#' 
#' @param verbose (logical) print messages
#' @examples getExamplePathways()
#' @return (char) Path to local cached copy of GMT file
#' or initial download is required 
#' @export
getExamplePathways <- function(verbose = FALSE) {
    pathwayURL <- paste("http://download.baderlab.org/EM_Genesets/", "January_24_2016/Human/symbol/", 
        "Human_AllPathways_January_24_2016_symbol.gmt", sep = "")
    bfc <- .get_cache()
    rid <- bfcquery(bfc, "Example_Pathways", "rname")$rid
    if (!length(rid)) {
        if (verbose) 
            message(paste("Downloading example pathway definitions (only ", "required once)", 
                sep = ""))
        rid <- names(bfcadd(bfc, "Example_Pathways", pathwayURL))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid))) 
        bfcdownload(bfc, rid)
    
    bfcrpath(bfc, rids = rid)
}
