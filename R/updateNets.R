#' Synchronize patient set in sample table and network table.
#'
#' @details This function is useful in applications with highly missing
#' data or where each patient contributes data points not present in the
#' others; e.g. networks based on individual
#' patient CNVs, which are highly sparse. In such a scenario, any kind of
#' patient subsetting - for example, limiting to training samples - changes
#' the population of eligible networks for analysis. Networks that no longer
#' have samples, or that have one patient with the neighbour removed, have
#' to be excluded. This function updates networks and patients so that 
#' each network contains at least two patients and only patients in 
#' networks are retained. In other words, it keeps pheno_DF and p_net in 
#' sync.
#' @param p_net (matrix) rows are patients, columns are networks.
#' a[i,j] = 1 if patient i occurs in network j, else 0.
#' @param pheno_DF (data.frame) patient ID and STATUS. 
#' @param writeNewNets (logical) if TRUE writes new networks to 
#' \code{newNetDir}.
#' @param oldNetDir (char) path to directory with networks to be updated
#' @param newNetDir (char) path to directory where updated networks are
#' to be written
#' @param verbose (logical) print messages
#' @param ... passed to pruneNets()
#' @return list with updated p_net and pheno_DF. pheno_DF will contain IDs
#' in the updated p_net. p_net will contain only those networks with 
#' 2+ patients and those patients present in 1+ network.
#' @export
#' @examples
#' data(npheno)
#' netDir <- sprintf('%s/extdata/example_nets',path.package('netDx'))
#' netmat <- countPatientsInNet(netDir,dir(netDir,pattern='txt$'), npheno[,1])
#' x <- updateNets(netmat, npheno,writeNewNets=FALSE)
updateNets <- function(p_net, pheno_DF, writeNewNets = TRUE, oldNetDir, 
		newNetDir, verbose = TRUE, ...) {
    idx <- which(colSums(p_net) >= 2)
    p_net <- p_net[, idx]
    idx <- which(rowSums(p_net) >= 1)
    p_net <- p_net[idx, ]
    if (verbose) {
        message("Update: (num patients) x (num networks)")
        print(dim(p_net))
    }
    
    # training samples are only those that occur in label-enriched networks
    pheno_DF <- pheno_DF[which(pheno_DF$ID %in% rownames(p_net)), ]
    
    if (writeNewNets) {
        pruneNets(oldNetDir, newNetDir, filterNets = colnames(p_net), 
					filterIDs = rownames(p_net),  ...)
    }
    
    return(list(p_net = p_net, pheno_DF = pheno_DF))
}
