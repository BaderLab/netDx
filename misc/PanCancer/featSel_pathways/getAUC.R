# gets AUCROC and AUCPR from a table
getAUC <- function(dat,delta=0.5) {
aucroc <- NA; aucpr <- NA
if (nrow(dat)>0) {
		sc <-(dat$SURVIVEYES_SCORE-dat$SURVIVENO_SCORE) + delta
		sc <- (sc-(-1))/2 # normalize from 0 to 1
    pred <- prediction(sc,dat$STATUS=="SURVIVEYES")

###    c1 <- "SURVIVEYES" #numc[1]
###    tp <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS == c1)
###    tn <- sum(dat$STATUS==dat$PRED_CLASS & dat$STATUS != c1)
###    fp <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS != c1)
###    fn <- sum(dat$STATUS!=dat$PRED_CLASS & dat$STATUS == c1)

    # dummy score column needed for perfCalc()
    #tmp <- data.frame(score=0,tp=tp,tn=tn,fp=fp,fn=fn)
    auroc <- performance(pred, "auc")@y.values[[1]]
		x <- performance(pred,"tpr","fpr") 
		roc <- data.frame(tpr=x@x.values[[1]],x@y.values[[1]])
}
return(list(auroc=auroc,ROC=roc))
}
