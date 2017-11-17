#'Takes binary networks and appends "binary_" to as a prefix to the filename
#'
#' @param netDir (list) directory containing network files
#' @param netList (list) list of binary networks
#' @param pheno (data.frame) sample metadata, must have ID column
#' @examples

#' @export


updateBinaryNetlist <- function(netDir,netList,pheno) {
  p <- countPatientsInNet(netDir,netList, pheno$ID)
  tmp <- updateNets(p,pheno,writeNewNets=FALSE)
  p <- tmp[[1]]
  pheno <- tmp[[2]]
  ret_vals <- list(p = p, pheno = pheno)
  return(ret_vals)
  }
