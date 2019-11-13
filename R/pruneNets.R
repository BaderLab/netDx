#' Prune interaction networks to keep only the networks and patients 
#' requested
#'
#' @details This function is crucial for patient data that is highly 
#' sparse; examples include patient CNVs indels, as opposed to full matrix
#' measures (gene expression, questionnaire data). Each step where the pool
#' of patients is subset - e.g. limiting feature selection only to patients
#' in training set - changes the set of networks that are eligible. 
#' Some networks may only contain test patients, while others may contain
#' a single edge between a training and a test patient. Upon subsetting,
#' such networks are no longer eligible for downstream use, such as 
#' feature selection. This function rewrites those subnetworks of the 
#' original networks that consist of eligible patients. 
#' @param oldDir (char) path to directory with original networks
#' @param newDir (char) path to output directory for pruned networks
#' @param filterNets (char) vector of networks to include. These should 
#' match filenames in \code{netDir}. Value of "*" results in pruning all
#' networks
#' @param filterIDs (char) patients to include in pruned networks. These
#' should match nodes in the input interaction networks
#' @param netSfx (char) suffix for network file names. Only used if 
#' \code{filterNets="*"}.
#' @param verbose (logical) print messages
#' @return (no value). Side effect of writing pruned network files to 
#' \code{newDir}
#' @examples
#' data(npheno)
#' netDir <- sprintf("%s/extdata/example_nets",path.package("netDx"))
#' pruneNets(netDir,"~/tmp",filterIDs=npheno[seq_len(10),],
#' 	netSfx="txt$")
#' @export
pruneNets <- function(oldDir,newDir,filterNets="*",filterIDs="*",
	netSfx="_cont.txt$",verbose=TRUE) {
	if (length(filterNets)==1) {
		if (filterNets=="*") {
			if (verbose) message("* Including all networks\n")
			fList <- dir(path=oldDir,pattern=netSfx)
			filterNets <- fList
		}
	} 
	if (verbose) message(sprintf("Limiting to %i networks\n", 
			length(filterNets)))

	if (!file.exists(newDir)) dir.create(newDir)

	if (length(filterIDs)==1) {
		if (filterIDs=="*") {	# keep all patients
			message("* Including all patients\n")
			for (f in filterNets) {
				oldf <- sprintf("%s/%s",oldDir,f)
				newf <- sprintf("%s/%s", newDir,f)
				file.copy(oldf,newf)
			}
		}
	} else {
		if (verbose) message(sprintf("Limiting to %i patients\n", 
					length(filterIDs)))
		for (f in filterNets) {
			dat <- read.delim(sprintf("%s/%s",oldDir,f),
							  sep="\t",header=FALSE,as.is=TRUE)

			# both nodes of edge should be eligible
			idx <- intersect(which(dat[,1]%in% filterIDs),
							 which(dat[,2]%in% filterIDs))

			write.table(dat[idx,],file=sprintf("%s/%s",newDir,f),
						sep="\t",col=FALSE,row=FALSE,quote=FALSE)
		}
	}
}
