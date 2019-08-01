# functions to download and update GeneMANIA jar used for network integration

#' wrapper function for getting BiocFileCache associated with netDx package
#' 
#' @return BiocFileCache object associated with netDx
#' @import BiocFileCache
#' @importFrom rappdirs user_cache_dir
.get_cache <- function() {
    cache <- rappdirs::user_cache_dir(appname="netDx")
    BiocFileCache::BiocFileCache(cache)
}

#' download and update GeneMANIA jar file 
#' 
#' @param verbose (logical) print messages
#' @examples getGMjar_path()
#' @return (char) Path to local cached copy of GeneMANIA jar file..
#' or initial download is required 
#' @export
getGMjar_path <- function( verbose = FALSE ) {
    fileURL <- "http://download.baderlab.org/netDx/genemania-cytoscape-plugin-3.5.0.jar"

    bfc <- .get_cache()
    rid <- bfcquery(bfc, "GM_jar", "rname")$rid
    if (!length(rid)) {
     if( verbose )
		 message( "Downloading GeneMANIA jar file (only required once)" )
     	rid <- names(bfcadd(bfc, "GM_jar", fileURL))
    }
    if (!isFALSE(bfcneedsupdate(bfc, rid)))
    bfcdownload(bfc, rid)

    bfcrpath(bfc, rids = rid)
}
