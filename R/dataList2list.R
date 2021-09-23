#' Convert MultiAssayExperiment object to list and data.frame
#' 
#' @details Used by internal routines in netDx
#' @param dat (MultiAssayExperiment) Patient data and metadata
#' @param groupList (list) variable groupings used for feature construction. See groupList arg in buildPredictor().
#' @return (list) Keys are:
#' 1) assays: list of matrices, each corresponding to data from a particular
#' layer
#' 2) pheno: (data.frame) sample metadata 
#' @examples
#' data(xpr,pheno)
#' require(MultiAssayExperiment)
#' objlist <- list("RNA"=SummarizedExperiment(xpr))
#' mae <- MultiAssayExperiment(objlist,pheno)
#' groupList <- list(RNA=rownames(xpr))
#' dl <- dataList2List(mae,groupList)
#' summary(dl) 
#' @export
dataList2List <- function(dat,groupList) {

# convert assays to list of matrices, replacing assay-specific sample
# name with patient ID
exprs <- experiments(dat)
datList2 <- list()
for (k in seq_len(length(exprs))) {
	tmp <- exprs[[k]]
	df <- sampleMap(dat)[
			which(sampleMap(dat)$assay==names(exprs)[k]),]

	colnames(tmp) <- df$primary[match(df$colname,colnames(tmp))]
	if ("SimpleList" %in% class(tmp)){
		tmp <- as.matrix(assays(tmp)[[1]]) # convert to matrix
	} else if ("SummarizedExperiment" %in% class(tmp)){
		tmp <- as.matrix(assays(tmp)[[1]])
	}
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
