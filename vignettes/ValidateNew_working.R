suppressWarnings(suppressMessages(require(netDx)))
suppressMessages(library(curatedTCGAData))

# brca <- suppressMessages(curatedTCGAData("BRCA",c("mRNAArray"),FALSE, version="1.1.38"))
# staget <- sub("[abcd]","",sub("t","",colData(brca)$pathology_T_stage))
# staget <- suppressWarnings(as.integer(staget))
# colData(brca)$STAGE <- staget

# pam50 <- colData(brca)$PAM50.mRNA
# pam50[which(!pam50 %in% "Luminal A")] <- "notLumA"
# pam50[which(pam50 %in% "Luminal A")] <- "LumA"
# colData(brca)$pam_mod <- pam50

# tmp <- colData(brca)$PAM50.mRNA
# idx <- union(which(tmp %in% c("Normal-like","Luminal B","HER2-enriched")),
#              		which(is.na(staget)))
# pID <- colData(brca)$patientID
# tokeep <- setdiff(pID, pID[idx])
# brca <- brca[,tokeep,]

# # remove duplicate assays mapped to the same sample
# smp <- sampleMap(brca)
# samps <- smp[which(smp$assay=="BRCA_mRNAArray-20160128"),]
# notdup <- samps[which(!duplicated(samps$primary)),"colname"]
# brca[[1]] <- suppressMessages(brca[[1]][,notdup])

# message("* Holding out 10 samples per class for validation")
# set.seed(123)
# pam50 <- colData(brca)$pam_mod
# pheno <- colData(brca)
# idx_holdout <- c(sample(which(pam50 == "LumA"),10,F),
# 		 sample(which(pam50 == "notLumA"),10,F)
# 	)
# holdout <- brca[,rownames(pheno)[idx_holdout]]
# colData(holdout)$ID <- as.character(colData(holdout)$patientID)
# brca <- brca[,rownames(pheno)[setdiff(1:length(pam50),idx_holdout)]]

# pID <- as.character(colData(brca)$patientID)
# colData(brca)$ID <- pID
# colData(brca)$STATUS <- colData(brca)$pam_mod

# summary(brca)

# groupList <- list()

# # genes in mRNA data are grouped by pathways
# pathList <- readPathways(fetchPathwayDefinitions("January",2018))

# groupList[["BRCA_mRNAArray-20160128"]] <- pathList[1:3]
# # clinical data is not grouped; each variable is its own feature
# groupList[["clinical"]] <- list(
#       age="patient.age_at_initial_pathologic_diagnosis",
# 	   stage="STAGE"
# )

## -----------------------------------------------------------------------------
makeNets <- function(dataList, groupList, netDir, ...) {
  netList <- c() # initialize before is.null() check
  # make RNA nets (NOTE: the check for is.null() is important!)
  # (Pearson correlation)
  if (!is.null(groupList[["BRCA_mRNAArray-20160128"]])) {
    netList <- makePSN_NamedMatrix(dataList[["BRCA_mRNAArray-20160128"]],
        rownames(dataList[["BRCA_mRNAArray-20160128"]]),
          groupList[["BRCA_mRNAArray-20160128"]],
        netDir, verbose = FALSE,
          writeProfiles = TRUE, ...)
  }

  # make clinical nets (normalized difference)
  netList2 <- c()
  if (!is.null(groupList[["clinical"]])) {
    netList2 <- makePSN_NamedMatrix(dataList$clinical,
    rownames(dataList$clinical),
    groupList[["clinical"]], netDir,
    simMetric = "custom", customFunc = normDiff, # custom function
    writeProfiles = FALSE,
    sparsify = TRUE, verbose = TRUE, ...)
  }
  netList <- c(unlist(netList), unlist(netList2))
  return(netList)
}

# message("About to build predictor...")
# Sys.sleep(2)

# set.seed(42) # make results reproducible
# outDir <- paste(tempdir(),randAlphanumString(),
# 	"pred_output",sep=getFileSep())
# # set keepAllData=TRUE to not delete at the end of the predictor run.
# # This can be useful for debugging.
# out <- buildPredictor(
#       dataList=brca,groupList=groupList,
#       makeNetFunc=makeNets,
#       outDir=outDir, ## netDx requires absolute path
#       numSplits=2L,featScoreMax=2L,
#       featSelCutoff=1L,
#       numCores=1L,debugMode=FALSE,
#       logging="none")


# numSplits <- 2
# st <- unique(colData(brca)$STATUS)
# acc <- c()         # accuracy
# predList <- list() # prediction tables

# featScores <- list() # feature scores per class
# for (cur in unique(st)) featScores[[cur]] <- list()

# for (k in 1:numSplits) { 
# 	pred <- out[[sprintf("Split%i",k)]][["predictions"]];
# 	# predictions table
# 	tmp <- pred[,c("ID","STATUS","TT_STATUS","PRED_CLASS",
# 	                 sprintf("%s_SCORE",st))]
# 	predList[[k]] <- tmp 
# 	# accuracy
# 	acc <- c(acc, sum(tmp$PRED==tmp$STATUS)/nrow(tmp))
# 	# feature scores
# 	for (cur in unique(st)) {
# 	   tmp <- out[[sprintf("Split%i",k)]][["featureScores"]][[cur]]
# 	   colnames(tmp) <- c("PATHWAY_NAME","SCORE")
# 	   featScores[[cur]][[sprintf("Split%i",k)]] <- tmp
# 	}
# }

# ## ----eval=TRUE,fig.width=8,fig.height=10--------------------------------------
# predPerf <- plotPerf(predList, predClasses=st)

# featScores2 <- lapply(featScores, getNetConsensus)
# summary(featScores2)
# head(featScores2[["LumA"]])

# featSelNet <- lapply(featScores2, function(x) {
#     callFeatSel(x, fsCutoff=1, fsPctPass=0)
# })
# print(head(featScores2[["LumA"]]))

# blah <- lapply(featSelNet, function(x) {
# 	x <- gsub(".profile$","",x)
# 	x <- gsub("_cont.txt","",x)
# })

# message("* Classify held-out set using selected features")

# save(brca,holdout,groupList,featSelNet,outDir,file="brca_save.rda")

load("brca_save.rda")
colData(holdout)$STATUS <- colData(holdout)$pam_mod

outDir <- paste(tempdir(), randAlphanumString(), sep = getFileSep())
if (file.exists(outDir)) unlink(outDir, recursive = TRUE)
dir.create(outDir)

names(groupList[["BRCA_mRNAArray-20160128"]]) <- gsub(",", ".",
  names(groupList[["BRCA_mRNAArray-20160128"]]))
out <- predict(brca, holdout, groupList, featSelNet, makeNets,
  outDir, verbose = FALSE)

perf <- getPerformance(out, c("LumA", "notLumA"))
message(sprintf("AUROC=%1.2f", perf$auroc))
message(sprintf("AUPR=%1.2f", perf$aupr))
message(sprintf("Accuracy=%1.1f%%", perf$acc))

plotDir <- "/.mounts/labs/pailab/spai/sftp/share/netDx_dev"
require(plotrix)

pdf(sprintf("perf.pdf"))
plotPerf_multi(list(perf$rocCurve),
  plotTitle = sprintf("BRCA Validation: %i samples", nrow(colData(holdout))))
plotPerf_multi(list(perf$prCurve), plotType = "PR",
  plotTitle = sprintf("BRCA Validation: %i samples", nrow(colData(holdout))))
dev.off()

pdf("colormat.pdf")
tbl <- as.matrix(table(out[,c("STATUS","PRED_CLASS")]))
color2D.matplot((tbl/nrow(out))*100,show.values=TRUE,axes=FALSE,
	xlab="Predicted label", ylab="Actual label",par(las=1))
axis(1,at=seq_len(ncol(tbl))-0.5,labels=colnames(tbl))
axis(2,at=seq_len(nrow(tbl))-0.5,labels=rownames(tbl))
dev.off()

