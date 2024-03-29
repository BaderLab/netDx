#' predict patient labels 
#' 
#' @description Once a model is trained, this function is used to classify new patients using selected features
#' @param trainMAE (MultiAssayExperiment) patient data for training samples. Same as provided to buildPredictor()
#' @param testMAE (MultiAssayExperiment) new patient dataset for testing model. Assays must be the same as for trainMAE.
#' @param groupList (list) list of features used to train the model. Keys are data types, and values are lists for groupings within those datatypes.
#' e.g. keys could include {'clinical','rna','methylation'}, and values within 'rna' could include pathway names {'cell cycle', 'DNA repair'}, etc.,
#' featSel will be used to subset 
#' @param featSel (list) selected features to be used in the predictive model. 
#' keys are patient labels (e.g. "responder/nonresponder"), and values are feature names 
#' identified by running buildPredictor(). Feature names must correspond to names of groupList, from which they will be subset.
#' @param makeNetFunc (function) function to create PSN features from patient data. See makeNetFunc in buildPredictor() for details
#' @param impute (logical) if TRUE imputes train and test samples separately before creating features. Currently unsupported.
#' @param outDir (char) directory for results
#' @param verbose (logical) print messages
#' @param numCores (integer) number of CPU cores for parallel processing
#' @param JavaMemory (integer) memory in (Gb) used for each fold of CV
#' @param debugMode (logical) Set to TRUE for detailed messages. Used for debugging.
#' @return (data.frame) predicted patient similarities and labels
#' columns are:  1) ID, 2) STATUS (ground truth), 3) <label>_SCORE: similarity score for the corresponding label,
#' 4) PRED_CLASS: predicted class
#' @export
predict <- function(trainMAE, testMAE, groupList, featSel, makeNetFunc, outDir,
    impute = FALSE, verbose = FALSE, numCores = 1L, JavaMemory = 4L, debugMode = FALSE) {

  # input checks
  if (missing(trainMAE)) stop("trainMAE must be supplied.\n")
  if (missing(testMAE)) stop("testMAE must be supplied.\n")
  if (missing(groupList)) stop("groupList must be supplied.\n")
  if (length(groupList) < 1) stop("groupList must be of length 1+\n")
  if (class(featSel) != "list") stop("featSel must be a list with patient labels as keys, and selected features as values")
  if (missing(outDir)) stop("outDir must be supplied.\n")
  if (impute) stop("impute=TRUE is not supported in the current version of netDx. This will be implemented based on future user interest. Please contact Shraddha Pai <shraddha.pai@utoronto.ca> if this feature is required.")

  nm1 <- setdiff(names(groupList), "clinical")
  if (!is(trainMAE, "MultiAssayExperiment"))
    stop("trainMAE must be a MultiAssayExperiment")
  if (!is(testMAE, "MultiAssayExperiment"))
    stop("testMAE must be a MultiAssayExperiment")

  tmp <- unlist(lapply(groupList, class))
  not_list <- sum(tmp == "list") < length(tmp)

  names_nomatch <- any(!nm1 %in% names(trainMAE))
  if (!is(groupList, "list") || not_list || names_nomatch) {
    msg <- c("groupList must be a list of lists.",
  " Names must match those in trainMAE, and each entry should be a list",
  " of networks for this group.")
    stop(paste(msg, sep = ""))
  }

  for (nm in names(featSel)) {
    featSel[[nm]] <- sub("_cont.txt", "", sub(".profile", "", featSel[[nm]]))
  }
  # clean features
  fs <- unlist(featSel);
  names(fs) <- NULL
  gl <- c()
  for (k in names(groupList)) {
    m <- groupList[[k]] # dataset level
    gl <- c(gl, names(m)) # add entries within dataset level
  }

  if (sum(!fs %in% gl) > 0) {
    stop("One or more entry in featSelNet not found in groupList.")
  }

  # merging train-test for joint db
  trainList <- dataList2List(trainMAE,groupList)
  testList <- dataList2List(testMAE,groupList)

  ph <- trainList$pheno[, c("ID", "STATUS")]
  ph2 <- testList$pheno[, c("ID", "STATUS")]
  ph$TT_STATUS <- "TRAIN"
  ph2$TT_STATUS <- "TEST"

  message("* Merging metadata tables...")
  tryCatch({
    pheno <- rbind(ph, ph2)
  }, error = function(ex) {
    stop(paste("couldn't combine train and test pheno.",
            "check that they have identical columns in same order", sep = ""))
  })
  print(table(pheno[, c("STATUS", "TT_STATUS")]))


  message("* Merging assays ...")
  assays <- list()
  for (nm in names(trainList$assays)) {
    message(sprintf("\t%s", nm))
    tryCatch({
      assays[[nm]] <- cbind(trainList$assays[[nm]], testList$assays[[nm]])
    }, error = function(ex) {
      stop(sprintf(paste("Error while combining data type %s for train and test ",
            "samples. Have you checked that measures are identical for both?", sep = ""),
            nm));
    })
  }
  message("* Measuring similarity to each known class")
  subtypes <- unique(ph$STATUS)
  predRes <- list()
  for (g in subtypes) {
    if (verbose) message(sprintf("\t%s", g))
    pDir <- paste(outDir, g, sep = getFileSep())
    dir.create(pDir)

    netDir <- paste(pDir, "networks", sep = getFileSep())
    dir.create(netDir)
    pheno_id <- setupFeatureDB(pheno, netDir)

    if (verbose) message("Creating PSN")
    createPSN_MultiData(dataList = assays,
        groupList = groupList,
        pheno = pheno_id,
        netDir = netDir,
        customFunc = makeNetFunc,
        numCores = 1L,
        filterSet = featSel[[g]],
        verbose = verbose)
    dbDir <- compileFeatures(netDir,
      outDir = pDir,
      numCores = numCores,
      verbose = verbose,
      debugMode = debugMode)

    # run query for this class
    qSamps <- pheno$ID[which(pheno$STATUS %in% g & pheno$TT_STATUS %in% "TRAIN")]
    qFile <- paste(pDir, sprintf("%s_query", g), sep = getFileSep())

    message(sprintf("\t%s : %s training samples", g, prettyNum(length(qSamps), big.mark = ",")))
    writeQueryFile(qSamps, "all", nrow(pheno), qFile)
    if (verbose) message(sprintf("\t** %s: Compute similarity", g))
    resFile <- runQuery(dbDir$dbDir, qFile, resDir = pDir,
      JavaMemory = JavaMemory, numCores = numCores,
      verbose = verbose, debugMode = debugMode)
    predRes[[g]] <- getPatientRankings(sprintf("%s.PRANK", resFile), pheno, g)
  }

  predClass <- predictPatientLabels(predRes,
      verbose = verbose)
  out <- merge(x = pheno, y = predClass, by = "ID")
  acc <- sum(out$STATUS == out$PRED_CLASS) / nrow(out)
  message(sprintf("%s test patients", prettyNum(nrow(out), big.mark = ",")))
  message(sprintf("ACCURACY (N=%i test) = %2.1f%%",
      nrow(out), acc * 100))
  message("Confusion matrix")
  print(table(out[, c("STATUS", "PRED_CLASS")]))

  out <- out[, - which(colnames(out) == "TT_STATUS")]

  return(out)

}