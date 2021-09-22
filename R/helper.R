# CBW helper functions

#' Compiles performance and selected features for a trained model.
#'
#' @details This function is run after training a model using buildPredictor(). 
#' It takes patient input data, model output, and returns performance and selected features. 
#' @param dat (MultiAssayExperiment) input data
#' @param res (list) output of buildPredictor() function
#' @param featureSelCutoff (integer) cutoff score for feature selection.
#' A feature must have minimum of this score for specified fraction of splits 
#' (see featureSelPct) to pass.
#' @param featureSelPct (numeric between 0 and 1) cutoff percent for feature selection.
#' A feature must have minimum score of featureSelCutoff for featureSelPct of 
#' train/test splits to pass.
#' @returns list of results.
#' - selectedFeatures (list of character vectors): list, one per class
#' - performance (list of mixed datatypes) including mean accuracy (meanAccuracy), 
#' split-level accuracy (splitAccuracy), split-level AUROC (auroc),
#' split-level AUPR (splitAUR)
#' Side effect of plotting ROC curve if binary classifier
#' @export
getResults <- function(dat, res, featureSelCutoff=1L, 
    featureSelPct=0){

numSplits <- length(grep("^Split",names(res)))
st <- unique(colData(dat)$STATUS)
message(sprintf("Detected %i splits and %i classes", numSplits, length(st)))

acc <- c()         # accuracy
predList <- list() # prediction tables
featScores <- list() # feature scores per class
for (cur in unique(st)) featScores[[cur]] <- list()

# collect accuracy and feature scores
for (k in 1:numSplits) {
	pred <- res[[sprintf("Split%i",k)]][["predictions"]];
	# predictions table
	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
	                 sprintf("%s_SCORE",st))]
	predList[[k]] <- tmp
	# accuracy
	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
	# feature scores
	for (cur in unique(st)) {
	   tmp <- res[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
	}
}

# only plot ROC and PR curves 
auroc <- NULL; aupr <- NULL
if (length(st)==2) {
message("* Plotting performance")
predPerf <- plotPerf(predList, predClasses=st)
auroc <- unlist(lapply(predPerf, function(x) x$auroc))
aupr <- unlist(lapply(predPerf, function(x) x$aupr))
}

message("* Compiling feature scores and calling selected features")
feats <- callOverallSelectedFeatures(featScores, 
    featureSelCutoff = featureSelCutoff,
    featureSelPct = featureSelPct,
    cleanNames = TRUE
)

#### Enrichment map
###if (!is.null(pathwayList)){
###    message("* Pathway List detected - creating input for EnrichmentMap")
###    browser()
###}

return(list(
    selectedFeatures=feats$selectedFeatures,
    featureScores=feats$featScores,
    performance=list(meanAccuracy=mean(acc),
                    splitAccuracy=acc,
                    splitAUROC=auroc,
                    splitAUPR=aupr)
))

}

#' Wrapper to call selected features
#'
#' @details Calls features that are consistently high-scoring for predicting 
#' each class. The context for this is as follows: 
#' The original model runs feature selection over multiple splits of data
#' into train/test samples, and each such split generates scores for all features.
#' This function identifies features with scores that exceed a threshold for a fraction
#' of train/test splits; the threshold and fraction are both user-specified. This
#' function is called by the wrapper getResults(), which returns both the matrix of 
#' feature scores across splits and list of features that pass the user-specified cutoffs.
#' @param featScores (list of lists): matrix of feature scores across all splits, separated
#' by patient label. First level: patient labels. Second level: matrix of scores for 
#' corresponding label.
#' @param featureSelCutoff (integer) cutoff score for feature selection.
#' A feature must have minimum of this score for specified fraction of splits 
#' (see featureSelPct) to pass.
#' @param featureSelPct (numeric between 0 and 1) cutoff percent for feature selection.
#' A feature must have minimum score of featureSelCutoff for featureSelPct of 
#' train/test splits to pass.
#' @param cleanNames (logical) remove internal suffixes for human readability
#' @return (list) Feature scores for all splits, plus those passing selection for overall predictor
#' featScores: (matrix) feature scores for each split
#' selectedFeatures: (list) features passing selection for each class; one key per class
#' @export
callOverallSelectedFeatures <- function(featScores, featureSelCutoff, 
    featureSelPct, cleanNames=TRUE){
featScores2 <- lapply(featScores, getNetConsensus)
if (cleanNames) {
    featScores2 <- lapply(featScores2,function(x){
        x$PATHWAY_NAME <- sub(".profile","",x$PATHWAY_NAME)
        x$PATHWAY_NAME <- sub("_cont.txt","",x$PATHWAY_NAME)
        colnames(x)[1] <- "Feature"
        x
    })
}
featSelNet <- lapply(featScores2, function(x) {
    x <- callFeatSel(x, fsCutoff=featureSelCutoff, fsPctPass=featureSelPct)
})

return(list(
    featScores=featScores2,
    selectedFeatures=featSelNet
))
}

#' Wrapper to create input files for Enrichment Map
#'
#' @details An Enrichment Map is a network-based visualization of top-scoring pathway features
#' and themes. It is generated in Cytoscape. This script generates the input files needed
#' for Cytoscape to create an Enrichment Map visualization.
#' @param model (list) Output of training model, generated by running buildPredictor()
#' @param results (list) Model results. output of getResults()
#' @param pathwayList (list) output of readPathwayFile() used to make pathway-level feat ures for predictor
#' @param EMapMinScore (integer) minimum score for Enrichment Map
#' @param EMapMaxScore (integer) maximum score for Enrichment Map
#' @param EMapPctPass (numeric between 0 and 1) percent of splits for which feature must have score in range
#'  [EMapMinScore,EMapMaxScore] to be included for EnrichmentMap visualization
#' @param outDir (char) directory where files should be written
#' @return 
#' @export
makeInputForEnrichmentMap <- function(model,results,pathwayList,
    EMapMinScore=0L, EMapMaxScore=1L,
    EMapPctPass=0.5,outDir)
{
    featScores <- results$featureScores

message("* Creating input files for EnrichmentMap")
Emap_res <- getEMapInput_many(featScores,
    pathwayList,
    minScore=EMapMinScore,
    maxScore=EMapMaxScore,
    pctPass=EMapPctPass,
    model$inputNets,
    verbose=FALSE
)

gmtFiles <- list()
nodeAttrFiles <- list()

message("* Writing files for network visualization")
for (g in names(Emap_res)) {
    outFile <- paste(outDir,sprintf("%s_nodeAttrs.txt",g),sep=getFileSep())
    write.table(Emap_res[[g]][["nodeAttrs"]],file=outFile,
        sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
    nodeAttrFiles[[g]] <- outFile

    outFile <- paste(outDir,sprintf("%s.gmt",g),sep=getFileSep())
    conn <- suppressWarnings(
         suppressMessages(base::file(outFile,"w")))
    tmp <- Emap_res[[g]][["featureSets"]]
    gmtFiles[[g]] <- outFile

    for (cur in names(tmp)) {
        curr <- sprintf("%s\t%s\t%s", cur,cur,
            paste(tmp[[cur]],collapse="\t"))
        writeLines(curr,con=conn)
    }
close(conn)
}

return(list(GMTfiles=gmtFiles,NodeStyles=nodeAttrFiles))
}

#' get the integrated patient similarity network made of selected features
#' 
#' @details An integrated patient similarity network can be built using combined
#' top features for each patient class. Such a network is created by taking the union of selected features for
#' all patient labels, and aggregating pairwise edges for all of them using a user-specified function (aggFun).
#' The network is then pruned prior to visualization, using a user-specified fraction of strongest edges
#' (prune_pctX, prune_useTop). In addition, the user may quantify the distance between patients of the 
#' same class, relative to those of other classes, using Dijkstra distance (calcShortestPath flag).  
#' @param dat (MultiAssayExperiment) input data
#' @param groupList (list) feature groups, identical to groupList provided for buildPredictor()
#' @param makeNets (function) Function used to create patient similarity networks. Identical to 
#' makeNets provided to buildPredictor()
#' @param selectedFeatures (list) selected features for each class (key of list). This object is returned as
#' part of a call to getResults(), after running buildPredictor().
#' @param plotCytoscape (logical) If TRUE, plots network in Cytoscape.
#' Requires Cytoscape software to be installed and running on the computer
#' when the function call is being made.
#' @param aggFun (char) function to aggregate edges from different PSN (e.g. mean)
#' @param prune_pctX (numeric between 0 and 1) fraction of most/least 
#' edges to keep when pruning the integrated PSN for visualization.
#' Must be used in conjunction with useTop=TRUE/FALSE
#' e.g. Setting pctX=0.2 and useTop=TRUE will keep 20\% top edges
#' @param prune_useTop (logical) when pruning integrated PSN for visualization,
#' determines whether to keep strongest edges (useTop=TRUE) or weakest edges
#' (useTop=FALSE)
#' @param numCores (integer) number of cores for parallel processing
#' @param calcShortestPath (logical) if TRUE, computes weighted shortest path
#' Unless you plan to analyse these separately from looking at the shortest 
#' path violin plots or integrated PSN in Cytoscape, probably good to set to 
#' FALSE.
#' @return (list) information about the integrated network
#  1) patientSimNetwork_unpruned (matrix) full integrated 
#' similarity network
#' 2) patientDistNetwork_pruned (matrix) the network plotted in
#' Cytoscape. Also note that this is a dissimilarity network, 
#' so that more similar nodes have smaller edge weights
#' 3) colLegend (data.frame): legend for the patient network
#' plotted in Cytoscape. Columns are node labels (STATUS) and
#' colours (colour)
#' 6) outDir (char) value of outDir parameter
#' @export
getPSN <- function(dat, groupList, makeNets, selectedFeatures, plotCytoscape=FALSE,
    aggFun="MEAN", prune_pctX=0.30, prune_useTop=TRUE,numCores=1L,calcShortestPath=FALSE
    ){
topPath <- gsub(".profile","", unique(unlist(selectedFeatures)))
topPath <- gsub("_cont.txt","",topPath)

## create groupList limited to top features
g2 <- list();
for (nm in names(groupList)) {
	cur <- groupList[[nm]]
	idx <- which(names(cur) %in% topPath)
	message(sprintf("%s: %i features", nm, length(idx)))
	if (length(idx)>0) g2[[nm]] <- cur[idx]
}

message("* Making integrated PSN")
psn <- 
   plotIntegratedPatientNetwork(dat,
  groupList=g2, makeNetFunc=makeNets,
  aggFun=aggFun,
  prune_pctX=prune_pctX,
  prune_useTop=prune_useTop,
  numCores=numCores,
  calcShortestPath=calcShortestPath,
  showStats=FALSE,
  verbose=TRUE, 
  plotCytoscape=plotCytoscape)

return(psn)
}

#' Make confusion matrix
#' 
#' @details Creates a confusion matrix, a square matrix which indicates the fraction of times
#' patients in a class are correctly classified, versus misclassified as each of the other classes.
#' Here, the confusion matrix is computed once per train-test split and the average is displayed. 
#' For this reason, the fractions may not cleanly add up to 100%.
#' @param model (list) output of buildPredictor()
#' @return (list) confusion matrix for all train/test splits and final averaged matrix
#' Side effect of plotting the averaged matrix.
#' @importFrom plotrix color2D.matplot
#' @export
confusionMatrix <- function(model) {
    nmList <- names(model)[grep("Split",names(model))]
    cl <- sort(unique(model$Split1$STATUS))
    conf <- list()
    mega <- NULL
    for (nm in nmList){
        pred <- model[[nm]][["predictions"]][,c("ID","STATUS","TT_STATUS","PRED_CLASS")]
        m <- as.matrix(table(pred[,c("STATUS","PRED_CLASS")]))
        conf[[nm]] <- m/colSums(m)
        if (is.null(mega)) mega <- conf[[nm]] else mega <- mega + conf[[nm]]
    }
    
        mega <- mega / length(conf) # average
        mega <- round(mega*100,2)
        mega <- t(mega)
        metric <- "%% Accuracy"
        
        tbl <- table(model$Split1$predictions$STATUS)
        nm <- names(tbl); val <- as.integer(tbl)
        ttl <- sprintf("%s\n(N=%i)",rownames(mega),val[match(rownames(mega),nm)])

    par(mar=c(4,8,2,2))
    color2D.matplot(mega,show.values=TRUE, border="white", 
        #cs1=c(1,1,1),cs2=c(1,0.5,0),cs3=c(0,0.5,0), 
        extremes=c(1,2),
        axes=FALSE,        
        xlab="Predicted class",ylab="")
    axis(1,at=seq_len(ncol(mega))-0.5,labels=colnames(mega))
    axis(2,at=seq_len(ncol(mega))-0.5,labels=rev(ttl),las=2)
    title(sprintf("Confusion matrix: Accuracy (avg of %i splits)",length(conf)))

    return(list(splitWiseConfMatrix=conf, average=mega))
}

#' Plot tSNE
#' 
#' @details Plots tSNE of integrated patient similarity network using Rtsne
#' @param psn (matrix) Patient similarity network represented as adjacency
#' matrix (symmetric). Row and column names are patient IDs. Note that NA
#' values will be replaced by very small number (effectively zero).
#' @param pheno (data.frame) Patient labels. ID column is patient ID and 
#' STATUS is patient label of interest. tSNE will colour-code nodes by 
#' patient label.
#' @param ... Parameters for Rtsne() function.
#' @return (Rtsne) output of Rtsne call. Side effect of tSNE plot
#' @import ggplot2
#' @import Rtsne Rtsne
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggplot2 ggplot
#' @export
tSNEPlotter <- function(psn,pheno,...) {

message("* Making symmetric matrix")
symmForm <- suppressMessages(makeSymmetric(psn))
symmForm[which(is.na(symmForm))] <- .Machine$double.eps
message("* Running tSNE")
x <- Rtsne(symmForm,...)
dat <- x$Y
samps <- rownames(symmForm)
idx <- match(samps, pheno$ID)
if (all.equal(pheno$ID[idx],samps)!=TRUE) {
	stop("pheno IDs not matching psn rownames")
}
st <- pheno$STATUS[idx]

# to eliminate the "no visible binding for global variable" problem
y <- status <- NULL

message("* Plotting")
colnames(dat) <- c("x","y")
dat <- as.data.frame(dat,stringsAsFactors=TRUE)
dat$status <- as.factor(st)
p <- ggplot2::ggplot(dat,aes(x,y)) + geom_point(aes(colour=status))
p <- p + xlab("") + ylab("") + ggtitle("Integrated PSN - tSNE")
print(p)

return(x)
}