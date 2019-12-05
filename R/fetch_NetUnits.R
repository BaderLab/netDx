#' Find the units within features, that patients individually map to
#'
#' @details This function is useful for binary networks where patients
#' contributes to units of a set. For example, CNV networks where patient
#' overlaps with single features - here, genes - result in membership in a
#' set - here, pathways. 
#' @param pat_GR (GRanges) ranges that patients map to. The metadata 
#' columns are required to have:
#' 1. ID (char): patient IDs
#' 2. LOCUS_NAMES (char): comma-separated vector of units that patients
#' contributed 
#' @param netList (list) keys are networks, values are vectors of 
#' unit features within those networks (e.g. pathways and genes)
#' @param whichNets (char) which networks to fetch data for ; e.g. 
#' could be set of feature-selected networks
#' @param trackMapping_detail (logical) if TRUE, returns the mapping
#' from individual patient events to units to features.
#' If FALSE, simply groups all units that contribute to 'whichNets',
#' without tracking individual mappings. 
#' The former is useful for understanding the contribution of individual
#' events to the overall feature selection.
#' @param verbose (logical) print messages
#' @return Depends on trackMapping_detail value.
#' if FALSE, (data.frame) ID: Patient ID ; UNIT: comma-separated vector of
#' unit feature names contributed by the patient
#' if TRUE, (GRanges) same as input pat_GR, but with a metadata column 
#' added.
#' New column is named unit2feature. Contains semicolon-delimited pairs of 
#' genes and feature-selected pathways in which genes are members. 
#' Format: geneX@pathway1;geneX@pathway2;geneY@pathway3. 
#' Here, the CNV overlaps geneX, which belongs in pathway1 and pathway2, 
#' and  also geneY, which belongs in pathway3. If a CNV does not overlap 
#' any genes in any feature-selected pathways, the value for that CNV is 
#' 'NA'.
#' @export
#' @examples
#' data(cnv_GR,pathway_GR,pathwayList)
#' x <- getRegionOL(cnv_GR,pathway_GR)
#' y <- fetch_NetUnits(x,pathwayList, names(pathwayList))
#' y <- fetch_NetUnits(x,pathwayList, names(pathwayList),
#' \ttrackMapping_detail=TRUE)
fetch_NetUnits <- function(pat_GR, netList, whichNets, trackMapping_detail = FALSE, 
    verbose = FALSE) {
    netg <- NULL
    for (n in whichNets) {
        netg <- c(netg, netList[[n]])
    }
    netg <- unique(unlist(netg))
    message(sprintf("Total %i subfeatures\n", length(netg)))
    
    patID <- unique(pat_GR$ID)
    patUnit <- rep(NA, length(patID))
    patNets <- rep(NA, length(patID))
    ctr <- 1
    
    if (trackMapping_detail) {
        message("tracking detailed mapping. indexing unit membership\n")
        unit_mat <- matrix(0, nrow = length(whichNets), ncol = length(netg))
        colnames(unit_mat) <- netg
        rownames(unit_mat) <- whichNets
        ctr <- 1
        for (netName in whichNets) {
            idx <- which(netg %in% netList[[netName]])
            if (verbose) 
                message(sprintf("%s\n", netName))
            unit_mat[ctr, idx] <- 1
            ctr <- ctr + 1
        }
        
        ### c() map unit-to-feature for each structural variant
        outcol <- c()
        for (k in seq_len(length(pat_GR))) {
            myg <- unlist(strsplit(pat_GR$LOCUS_NAMES[k], ","))
            
            # which genes are in fs pathways?
            gene_col <- which(colnames(unit_mat) %in% myg)
            # get the pathways in which this gene occurs
            has_gene <- rowSums(unit_mat[, gene_col, drop = FALSE])
            idx <- which(has_gene > 0)
            
            # create a semi-colon delimited pair of 'gene:pathway' for this cnv
            if (length(idx) > 0) {
                str <- c()
                for (i in idx) {
                  feat <- rownames(unit_mat)[i]
                  pull_gene <- intersect(myg, netList[[feat]])
                  str <- c(str, paste(pull_gene, feat, sep = "@"))
                }
                str <- paste(str, collapse = ";")
                outcol <- c(outcol, str)
            } else {
                outcol <- c(outcol, "NA")
            }
        }
        
        pat_GR$unit2feature <- outcol
        out <- pat_GR
        
    } else {
        for (uq in patID) {
            myidx <- which(pat_GR$ID == uq)
            myg <- unlist(strsplit(pat_GR$LOCUS_NAMES[myidx], ","))
            
            # genes that overlap my CNV and that are in
            patUnit[ctr] <- paste(intersect(myg, netg), collapse = ",")
            
            tmp <- c()
            for (n in whichNets) {
                m <- intersect(myg, netList[[n]])
                if (length(m) > 0) 
                  tmp <- c(tmp, n)
            }
            patNets[ctr] <- paste(tmp, collapse = ",")
            ctr <- ctr + 1
        }
        out <- data.frame(ID = patID, FEATURE = patNets, UNIT = patUnit)
    }
    
    return(out)
}
