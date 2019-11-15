#' Return feature selected nets based on given criteria
#'
#' @details given the output of genNetScores.R and criteria for defining
#' feature-selected (FS) nets, returns subset of nets that pass criteria.
#' Net must score <fsCutoff> for at least <fsPctPass> % of splits, to be
#' considered feature-selected.
#' @param netScores (matrix) matrix of net scores
#' @param fsCutoff (integer) net must score at least this much in a split to
#'  "pass" the threshold
#' @param fsPctPass (numeric 0 to 1) net must pass at least this percent of
#' splits to be considered feature-selected
#' @return (char) names of nets that pass feature-selection
#' @examples
#' data(featScores)
#' passed <- lapply(featScores, function(x) {
#'    callFeatSel(x,10,0.7) # score 10/10 in >=70% of trials
#' })
#' print(passed)
#' @export
callFeatSel <- function(netScores,fsCutoff, fsPctPass) {
  fs_nets <- c()
  for (index in seq_len(nrow(netScores))){
    cur_pathway <- netScores[index,]
    pass_thresh <- length(which(cur_pathway >= fsCutoff))
    percent_pass <- pass_thresh/length(cur_pathway)
    if(percent_pass >= fsPctPass){
      fs_nets <- c(fs_nets, netScores[,1][index])
    }
  }
  return(fs_nets)
}
