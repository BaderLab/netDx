#' cleaner sparsification routine - faster, matrix-based version
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
#' then the strongest edge for each missing patient is added to the net. Note 
#' that this condition results in the total number of edges potentially 
#' exceeding EDGE_MAX
#' @param verbose (logical) print detailed messages, useful for debugging
#' @return writes SIF content to text file (node1,node2,edge weight)
#' @import reshape2
#' @importFrom utils write.table
#' @examples
#' m <- matrix(runif(500*500),nrow=500)
#' y <- sparsify2(m)
#' @examples 
#' m <- matrix(runif(500*500),nrow=500)
#' y <- sparsify2(m)
#' @export
sparsify3 <- function(W, outFile = "tmp.txt", cutoff = 0.3, maxInt = 50, EDGE_MAX = Inf, 
    includeAllNodes = TRUE, verbose = TRUE) {
    
    if (maxInt > ncol(W)) 
        maxInt <- ncol(W)
    
    if (!is(W, "matrix")) 
        W <- as.matrix(W)
    W[which(is.na(W))] <- .Machine$double.eps  # don't allow missing values
    diag(W) <- NA
    mytop <- cbind(colnames(W), colnames(W)[apply(W, 1, which.max)], apply(W, 1, 
        max, na.rm = TRUE))
    # don't want same patient edge twice, nor self-similarity
    W[upper.tri(W, diag = TRUE)] <- NA
    W[W < cutoff] <- NA
    maxind <- min(ncol(W), maxInt)
    
    # effectively empty out the slots that are not the top interactions create a
    # 'switch off' matrix with NA in non-top edges
    W_order <- t(apply(W, 1, order, decreasing = TRUE, na.last = TRUE))
    W_order[which(W_order > maxInt)] <- NA
    W_order[which(W_order <= maxInt)] <- .Machine$double.eps
    W2 <- W + W_order  # NA for non-top edges, unchanged for top edges
    mmat <- na.omit(melt(W2, varnames = names(dimnames(W2))))
    
    maxEdge <- nrow(mmat)
    if (!is.infinite(EDGE_MAX)) {
        if (maxEdge > EDGE_MAX) 
            maxEdge <- EDGE_MAX
        mmat <- mmat[seq_len(maxEdge), ]
    }
    
    # we should guarantee an edge from all patients- in this case the edge_max would
    # be violated unless we come up with a better rule message('past this\n')
    if (includeAllNodes) {
        mmat[, 1] <- as.character(mmat[, 1])
        mmat[, 2] <- as.character(mmat[, 2])
        univ <- c(mmat[, 1], mmat[, 2])
        missing <- setdiff(rownames(W), univ)
        # message(sprintf('missing = { %s }\n',paste(missing, collapse=',')))
        if (length(missing) > 0) {
            message(sprintf(paste("Sparsify2: found %i missing patients; ", "adding strongest edge\n", 
                sep = ""), length(missing)))
            for (k in missing) {
                # add the strongest edge for the patient
                tmp <- mytop[which(mytop[, 1] %in% k), ]
                x <- as.numeric(tmp[3])
                if (x < cutoff) {
                  message("\tMissing edge is below cutoff; setting to cutoff\n")
                  x <- cutoff
                }
                mmat <- rbind(mmat, c(k, tmp[2], x))
            }
        }
    }
    
    head(mmat)
    mmat <- na.omit(mmat)  # boundary case where cutoff exceeds net max
    mmat[, 3] <- as.numeric(mmat[, 3])
    mmat[, 3] <- round(mmat[, 3], digits = 4)
    write.table(mmat, file = outFile, sep = "\t", col.names = FALSE, row.names = FALSE, 
        quote = FALSE)
    return(mmat)
}

