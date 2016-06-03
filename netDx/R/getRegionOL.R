#' Returns overlapping named ranges for input ranges
#'
#' @details Given a set of query GRanges, and a subject list-of-GRanges,
#' updates the query with a column 'LOCUS_NAMES' containing the names of
#' ranges overlapped by the query. One application is to map structural
#' variants, such as CNVs, to genes in pathways of interest. In this 
#' scenario \code{gr} would contain the patient CNVs, and \code{rngList}
#' would be a list of GenomicRanges objects, one per cellular pathway.
#' @param gr (GRanges) query ranges
#' @param rngList (list) keys are names, and values are GRanges, each range
#' of which has a name (in "name" column). Note: It is faster to provide
#' a list of length 1 ; if the list is long, combining into a single GRanges
#' object could prove slow.
#' @return (GRanges) query ranges with the added column "LOCUS_NAMES". 
#' Where a range overlaps with multiple loci, the names are reported as a 
#' comma-separated vector
#' @examples
#' data(TCGA_mini,pathway_GR)
#' x <- getRegionOL(cnv_GR,pathway_GR)
#' @export
getRegionOL <- function(gr, rngList) {
	rng <- GRanges()
	for (k in 1:length(rngList)) {
		cur <- rngList[[k]]
		seqlevels(rng) <- unique(c(seqlevels(rng), seqlevels(cur)))
		rng <- c(rng, cur)
	}

	rng <- rng[which(as.character(seqnames(rng)) %in% 
					 as.character(seqlevels(gr)))]
	seqlevels(rng) <- seqlevels(gr)

	ol <- as.matrix(findOverlaps(gr, rng))

	# could be made more efficient.
	ol_nm <- rng$name[ol[,2]]
	LOCUS_NAMES <- rep("", length(gr))
	for (k in unique(ol[,1])) {
		idx <- which(ol[,1] == k)
		LOCUS_NAMES[k] <- paste(unique(ol_nm[idx]),collapse=",")
	}
	gr$LOCUS_NAMES=LOCUS_NAMES

	gr
}
