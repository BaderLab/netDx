#' Counts number of (+,+) and (+,-) interactions in a set of networks
#' 
#' @param inFiles (char) path to interaction networks to process
#' @param plusID (char) IDs of + nodes
#' @param minusID (char) IDs of - nodes
#' @param tmpDir (char) path to dir where temporary files can be stored
#' @param enrType (char) see getEnr.R
#' @import bigmemory
#' @return (matrix) two columns, one row per network 
#' If \code{enrType="binary"}, number of (+,+) and other interactions
#' Otherwise if \code{enrType="corr"} mean edge weight of (+,+) edges and
#' of other edges
#' @export
countIntType_batch <- function(inFiles,plusID, minusID,tmpDir="/tmp",
	   enrType="binary"){
	bkFile <- sprintf("%s/tmp.bk",tmpDir)
	if (file.exists(bkFile)) file.remove(bkFile)
	out <- big.matrix(NA, nrow=length(inFiles),ncol=2,
					  type="double",backingfile="tmp.bk",
					  backingpath=tmpDir,
					  descriptorfile="tmp.desc")
	foreach (k=1:length(inFiles)) %dopar% {
		m <- bigmemory::attach.big.matrix(
				sprintf("%s/tmp.desc",tmpDir))
		if (enrType == "binary")
			m[k,]	<- countIntType(inFiles[k], plusID,minusID)
		else if (enrType == "corr")
			m[k,]	<- getCorrType(inFiles[k],plusID,minusID)
	}
	
	out	<- as.matrix(out)
	unlink(sprintf("%s/tmp.bk",tmpDir))
	unlink(sprintf("%s/tmp.desc",tmpDir))
	return(out)
}
