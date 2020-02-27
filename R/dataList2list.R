#' Convert MultiAssayExperiment object to list and data.frame
#' 
#' @details Used by internal routines in netDx
#' @param dat (MultiAssayExperiment) 
#' @return (list) Keys are:
#' 1) assays: list of matrices, each corresponding to data from a particular
#' layer
#' 2) pheno: (data.frame) sample metadata 
#' @export
dataList2List <- function(dat) {

# convert assays to list of matrices, replacing assay-specific sample
# name with patient ID
exprs <- experiments(dat)
datList2 <- list()
for (k in seq_len(length(exprs))) {
	tmp <- exprs[[k]]
	df <- sampleMap(dat)[
			which(sampleMap(dat)$assay==names(exprs)[k]),]

	colnames(tmp) <- df$primary[match(df$colname,colnames(tmp))]
	tmp <- as.matrix(assays(tmp)[[1]]) # convert to matrix
	datList2[[names(exprs)[k]]]<- tmp	
}
if ("clinical" %in% names(groupList)) {
	tmp <- colData(dat)
	vars <- unique(unlist(groupList[["clinical"]]))
	datList2[["clinical"]] <- t(as.matrix(tmp[,vars,drop=FALSE]))
}

pheno_all <- colData(dat)
pheno_all <- as.data.frame(pheno_all)

out <- list(
	assays=datList2,
	pheno=pheno_all)
}
