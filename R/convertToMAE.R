#' Wrapper that converts an input list into a MultiAssayExperiment object
#' 
#' @details This function takes in a list of key-value pairs (keys: data types,
#' values: matrices/dataframes) and calls the necessary functions from the
#' MultiAssayExperiment package to incorporate the values from the input list 
#' into a MultiAssayExperiment object, transforming the values according to the 
#' keys. If duplicate sample names are found in the assay data, only the first
#' instance is kept.
#' @param dataList  (list) input key-value pairs (keys: data types, values: 
#' data in the form of matrices/dataframes); must have a key-value pair that
#' corresponds to patient IDs/metadata labelled pheno.
#' @return MAE (MultiAssayExperiment) data from input list incorporated into a
#' MultiAssayExperiment object, compatible with further analysis using the 
#' netDx algorithm.
#' @examples
#' data(xpr, pheno)
#' 
#' # Generate random proteomic data
#' prot <- matrix(rnorm(100*20), ncol=20)
#' colnames(prot) <- sample(pheno$ID, 20)
#' rownames(prot) <- sprintf("protein%i",1:100)	
#' 
#' myList <- list(rna = xpr, proteomic = prot, pheno = pheno)
#' 
#' MAE <- convertToMAE(myList)
#' @export
convertToMAE <- function(dataList) {
  
  # Check input data:
  if (class(dataList) != "list") {
    stop("dataList must be a list. \n")
  }
  if (is.null(dataList$pheno)) {
    stop("dataList must have key-value pair labelled pheno.\n")
  }
  if (length(dataList) == 1) {
    stop("dataList must have assay data to incorporate into a 
         MultiAssayExperiment object")
  }
  
  # Note that a MultiAssayExperiment object requires an ExperimentList and 
  # colData (sampleMap optional if each assay uses the same colnames)
  
  # Possible elements for ExperimentList:
  # - base::matrix (gene expression, microRNA, metabolomics, microbiome data)
  # - SummarizedExperiment::SummarizedExperiment (same as matrix, but capable
  #   of storing additional assay-level metadata)
  # - Biobase::ExpressionSet (legacy representation, use SummarizedExperiment)
  # - SummarizedExperiment::RangedSummarizedExperiment (range-based datasets; 
  #   gene expression, methylation, data types that refer to genomic positions)
  # - RaggedExperiment::RaggedExperiment (range-based datasets; copy number and
  #   mutation data, measurements by genomic positions)
  
  # Assumes that pheno is a DataFrame (or coerceable to be a DataFrame)
  patientPheno <- dataList$pheno
  
  # Generate ExperimentList from input dataList
  tmp <- NULL
  track <- c()
  datType <- names(dataList)
  for (k in 1:length(dataList)) {
    # For key-value pairs that aren't labelled pheno, transform into 
    # objects compatible with input into MultiAssayExperiment object
    if (names(dataList[k]) != "pheno") {
      
      # Remove duplicated columns (we keep the first column) in the assay data
      if (sum(duplicated(colnames(dataList[[k]]))) != 0) {
        dataList[[k]] <- dataList[[k]][,!duplicated(colnames(dataList[[k]]))]
      }
      
      # Assumes that data is of matrix class 
      # *(maybe implement matrix conversion into SummarizedExperiment in future)
      track <- c(track, k)
      tmp <- c(tmp, list(dataList[[k]]))
    }
  }
  names(tmp) <- datType[track]
  
  MAE <- MultiAssayExperiment(experiments = tmp, colData = patientPheno)
  
  return(MAE)
}