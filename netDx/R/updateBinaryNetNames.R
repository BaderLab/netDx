#'Takes binary networks and appends "binary_" to as a prefix to the filename
#'
#' @param netDir (list) directory containing network files
#' @param netList (list) list of binary networks to have their filenames udpated
#' @examples

#' @export

updateBinaryNetNames <- function(netDir,netList) {
  updatedNetList <- c()
  for(curNet in netList){
    curNetPath <- sprintf('%s/%s', netDir, curNet)
    newNetName <- paste('binary_',curNet, sep='')
    newNetPath <- sprintf('%s/%s', netDir, newNetName)
    file.rename(curNetPath, newNetPath)
    updatedNetList <- c(updatedNetList, newNetName)
  }
  return(updatedNetList)
}
