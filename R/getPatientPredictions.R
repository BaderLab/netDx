#' Calculates patient-level classification accuracy across train/test splits
#'
#' @details Takes all the predictions across the different train/test splits,
#' and for each patient, generates a score indicating how many times they were
#' classified by netDx as belonging to each of the classes. The result is that
#' we get a measure of individual classification accuracy across the different
#' train/test splits.
#' @param predFiles (char) vector of paths to all test predictions
#' (e.g. 100 files for a 100 train/test split design).
#' Alternately, the user can also  provide a single directory name, and allow
#' the script to retrieve prediction files.
#' Format is "rootDir/rngX/predictionResults.txt"
#' @param pheno (data.frame) ID=patient ID, STATUS=ground truth (known class
#' label). This table is required to get the master list of all patients, as
#' not every patient is classified in every split.
#' @param plotAccuracy (logical) if TRUE, shows fraction of times
#' patient is misclassified, using a dot plot
#' @return (list) of length 2.
#' 1) (data.frame) rows are patients, (length(predFiles)+2) columns.
#' Columns 1:length(predFiles): Predicted labels for a given split (NA if
#' patient was training sample for the split).
#' Column (length(predFiles)+1):
#' split, value is NA. Columns are : ID, REAL_STATUS, predStatus1,...
#' predStatusN.
#' Side effect of plotting a dot plot of % accuracy. Each dot is a patient, and
#' the value is "% splits for which patient was classified correctly".
#' @examples
#' inDir <- sprintf("%s/extdata/KIRC_output",
#'    path.package("netDx.examples"))
#' require(netDx.examples)
#' data(KIRC_pheno)
#' all_rngs <- list.dirs(inDir, recursive = FALSE)
#' all_pred_files <- unlist(lapply(all_rngs, function(x) {
#'     paste(x, "predictionResults.txt", sep = "/")}))
#' pred_mat <- getPatientPredictions(all_pred_files, KIRC_pheno)
#' @import ggplot2
#' @export
getPatientPredictions <- function(predFiles,pheno,plotAccuracy=FALSE) {
  if(length(predFiles) == 1){
		cat("predFiles is of length 1. Assuming directory\n")
    all_rngs <- list.dirs(predFiles, recursive = FALSE)
		all_rngs <- all_rngs[grep("rng",all_rngs)]
    predFiles <- unlist(lapply(all_rngs, function(x) {
			paste(x, "predictionResults.txt", sep = "/")
		}))
  } else {
		cat("predFiles is of length > 1. Assuming filenames provided\n")
	}

  output_mat <- matrix(NA, nrow=nrow(pheno), ncol=length(predFiles)+2)

  patient_list <- list()
  for (cur_pat in pheno$ID) patient_list[[cur_pat]] <- c()


	uq_mat <- matrix(NA, nrow=nrow(pheno),ncol=length(predFiles))
	rownames(uq_mat) <- pheno$ID
	for (ctr in 1:length(predFiles)){
		curFile <- predFiles[ctr]
    dat 		<- read.delim(curFile,sep="\t",header=TRUE,as.is=TRUE)
		for (k in 1:nrow(dat)) {
			uq_mat[which(rownames(uq_mat)==dat$ID[k]),ctr] <- dat$PRED_CLASS[k]
		}
	}
uq_mat <- as.data.frame(uq_mat)
pctCorr <- c()
for (k in 1:nrow(uq_mat)) {
	testCt <- sum(!is.na(uq_mat[k,]))
	cur <- sum(uq_mat[k,] == pheno$STATUS[k],na.rm=TRUE)/testCt
	pctCorr <- c(pctCorr,cur*100)
}
uq_mat <- cbind(uq_mat, pheno$STATUS,pctCorr)

spos <- gregexpr("\\/",predFiles)
# get the name of the iteration (rngX) assuming directory structure
# rngX/pathway_CV_score.txt
fNames <- lapply(1:length(spos), function(x) {
  n <- length(spos[[x]])
	y <- substr(predFiles[x], spos[[x]][n-1]+1,spos[[x]][n]-1)
	y
})
fNames <- unlist(fNames)
out <- uq_mat;
rownames(out) <- pheno$ID
colnames(out) <- c(fNames,"STATUS","pctCorrect")

if (plotAccuracy) {
p <- ggplot(out,aes(x=pctCorrect))+ geom_dotplot()
p <- p + ggtitle(sprintf("Patient-level classification accuracy (N=%i)",
		length(predFiles)))
p <- p + theme(axis.text=element_text(size=13),
		axis.title=element_text(size=13))
print(p)
return(list(predictions=output_mat,plot=p))
} else 

return(list(predictions=output_mat))

}
