#' cleaner sparsification routine
#' 
#' @details Sparsifies similarity matrix to keep strongest edges.
#' Sets diagonal and edges < cutoff to NA. Keeps strongest maxInt edges
#' per node. Ties are ignored. Keeps a max of EDGE_MAX edges in the network.
#' @param W (matrix) similarity matrix
#' @param outFile (char) path to file to write sparsified network
#' @param cutoff (numeric) edges with weight smaller than this are set to NA
#' @param maxInt (numeric) max num edges per node. 
#' @param EDGE_MAX (numeric) max num edges in network
#' @param includeAllNodes (logical) if TRUE, ensures at least one edge is 
#' present for each patient. This feature is required when sparsification 
#' excludes test patients that are required to be classified. If the 
#' sparsification rules exclude all edges for a patient and this flag is set, 
#' then the strongest edge for each missing patient is added to the net. 
#' Note that this condition results in the total number of edges potentially 
#' exceeding EDGE_MAX
#' @param verbose (logical) print detailed messages, useful for debugging
#' @return writes SIF content to text file (node1,node2,edge weight)
#' @import reshape2
#' @export
sparsify2 <- function(W, outFile = "tmp.txt", cutoff = 0.3, maxInt = 50, EDGE_MAX = 1000, 
    includeAllNodes = TRUE, verbose = TRUE) {
    
    if (verbose) 
        message(sprintf(paste("sparsify2:maxInt=%i;EDGE_MAX=%1.2f;", "cutoff=%1.2e;includeAllNodes=%s", 
            sep = ""), maxInt, EDGE_MAX, cutoff, includeAllNodes))
    
    if (maxInt > ncol(W)) 
        maxInt <- ncol(W)
    
    # don't want same patient edge twice, nor self-similarity
    W[upper.tri(W, diag = TRUE)] <- NA
    W[W < cutoff] <- NA
    x <- list()
    for (i in seq_len(nrow(W))) {
        x[[i]] <- sort(W[i, ], decreasing = TRUE, na.last = TRUE)
    }
    # message('past sorting\n')
    names(x) <- rownames(W)
    for (k in seq_len(length(x))) {
        # print(k)
        cur <- x[[k]]
        tokeep <- names(cur)[seq_len(min(length(cur), maxInt))]
        W[k, which(!colnames(W) %in% tokeep)] <- NA
    }
    # message('got past b\n')
    mmat <- na.omit(melt(W))
    mmat <- mmat[order(mmat[, 3], decreasing = TRUE), ]
    
    if (!is.infinite(EDGE_MAX)) {
        maxEdge <- nrow(mmat)
        if (maxEdge > EDGE_MAX) 
            maxEdge <- EDGE_MAX
        mmat <- mmat[seq_len(maxEdge), ]
    }
    
    # we should guarantee an edge from all patients- in this case the edge_max would
    # be violated unless we come up with a better rule
    if (includeAllNodes) {
        mmat[, 1] <- as.character(mmat[, 1])
        mmat[, 2] <- as.character(mmat[, 2])
        univ <- c(mmat[, 1], mmat[, 2])
        missing <- setdiff(rownames(W), univ)
        # message(sprintf('missing = { %s }\n',paste(missing, collapse=',')))
        if (length(missing) > 0) {
            message(sprintf(paste("Sparsify2: ", "found %i missing patients; adding strongest edge\n", 
                sep = ""), length(missing)))
            for (k in missing) {
                # add the strongest edge for the patient
                tmp <- x[[k]]
                if (is.na(tmp[1])) {
                  message(paste("\tMissing edge is below cutoff; ", "setting to cutoff\n", 
                    sep = ""))
                  tmp[1] <- cutoff
                }
                mmat <- rbind(mmat, c(k, names(tmp)[1], tmp[1]))
            }
        }
    }
    head(mmat)
    mmat[, 3] <- as.numeric(mmat[, 3])
    mmat[, 3] <- round(mmat[, 3], digits = 4)
    write.table(mmat, file = outFile, sep = "\t", col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    return(mmat)
}

