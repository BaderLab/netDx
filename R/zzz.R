.onLoad <- function(libname, pkgname) {
  options(java.parameters = c("-Xmx10G"))
  message("Checking if Java runtime is installed  ...")
  x <- system2("java", "--version", stdout = TRUE, stderr = NULL)
  if (length(x)>0) {
    message("Java detected.")
    verNum <- as.integer(strsplit(strsplit(x[1], " ")[[1]][2], "\\.")[[1]][1])
    if (verNum > 16) {
      stop(paste("Incorrect Java version.\n",
        " Your Java version is ", verNum, ".",
        " netDx requires Java 16 or earlier to run.",
        " Please see https://github.com/RealPaiLab/netDx/blob/master/README.md for instructions on how to do this.",
        sep="", collapse=""))
    }
  } else {
      stop("Java not detected. Please install before proceeding.")
  }
}