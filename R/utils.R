
#' platform-specific file separator
#' 
#' @description Returns OS-specific file separator
#' @return (char) "\\" if Windows, else "/"
#' @export
getFileSep <- function(){
  if (.Platform$OS.type=="windows") return("\\")
  else return(.Platform$file.sep)
}
