#' Map named ranges to corresponding set of named ranges
#'
#' @details Example application is when we have named ranges each
#' corresponding to genes or regulatory elements, and we wish to group
#' these ranges based on metabolic pathway.
#' @param gr (GRanges) named ranged to be grouped
#' @param rangeList (list) sets of range names
#' @param verbose (logical) print detailed messages
#'
#' @return RangeList. keys are names of \code{rangeList}, values are GRanges
#' @import GenomicRanges
#' @examples 
#' data(genes_mini,pathwayList); 
#' gene_GR<- GRanges(genes$chrom,
#'	IRanges(genes$txStart,genes$txEnd),name=genes$name2)
#' path_GRList <- mapNamedRangesToSets(gene_GR,pathwayList)
#' @export
mapNamedRangesToSets <- function(gr,rangeList,verbose=FALSE){ 
    out <- list()
    for (nm in names(rangeList)) {
        my_gr <- gr[which(gr$name %in% rangeList[[nm]])]
        if (verbose) cat(sprintf("%s: %i ranges\n",nm,length(my_gr)))
        out[[nm]] <- my_gr
    }
    out
}
