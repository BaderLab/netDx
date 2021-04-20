#' performance metrics for model
#' @param res (data.frame) result from predicting labels on held-out test set. output of predict() function. 
#' columns include ID, STATUS (ground truth) and PRED_CLASS (predicted label)
#' @param predClasses (character) patient labels used by classifier
#' @return (list)
#' 1) rocCurve: ROCR performance object for ROC curve
#' 2) prCurve: ROCR performance object for PR curve
#' 3) auroc: Area under ROC curve
#' 4) aupr: Area under PR curve
#' 5) accuracy: Accuracy
#' @import ROCR
#' @export
getPerformance <- function(res, predClasses) {

  # given output of performance('precall') compute AUC-PR
  prauc <- function(res) {
    x <- res@x.values[[1]] # recall
    y <- res@y.values[[1]] # precision

    # remove NAN
    idx <- which(is.nan(y))
    if (any(idx)) {
      x <- x[-idx]
      y <- y[-idx]
    }

    pracma::trapz(x, y)
  }

  pred_col1 <- sprintf("%s_SCORE", predClasses[1])
  pred_col2 <- sprintf("%s_SCORE", predClasses[2])

  idx1 <- which(colnames(res) == pred_col1)
  idx2 <- which(colnames(res) == pred_col2)
  pred <- ROCR::prediction(res[, idx1] - res[, idx2],
            res$STATUS == predClasses[1])

  st <- res$STATUS
  c1 <- predClasses[1]
  tp <- sum(res$STATUS == res$PRED_CLASS & res$STATUS == c1)
  tn <- sum(res$STATUS == res$PRED_CLASS & res$STATUS != c1)
  fp <- sum(res$STATUS != res$PRED_CLASS & res$STATUS != c1)
  fn <- sum(res$STATUS != res$PRED_CLASS & res$STATUS == c1)
  # entire curves
  curRoc <- ROCR::performance(pred, "tpr", "fpr")
  curPr <- ROCR::performance(pred, "prec", "rec")
  tmp <- data.frame(score = 0, tp = tp, tn = tn, fp = fp, fn = fn)

  # statistic
  auroc <- ROCR::performance(pred, "auc")@y.values[[1]]
  aupr <- prauc(curPr)
  corr <- sum(res$STATUS == res$PRED_CLASS)
  acc <- (corr/nrow(res))*100

  return(list(rocCurve=curRoc,prCurve=curPr,auroc=auroc,aupr=aupr,accuracy=acc))
}