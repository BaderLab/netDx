# get iteration with performance closest to average 

datDir <- "/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_GBM/output/pruned_180204"

require(netDx)
predSet <- sprintf("%s/rng%i/all/cutoff9/predictionResults.txt",datDir,1:100)
perf <- plotPerf(predSet,predClasses=c("SURVIVEYES","SURVIVENO"))
auroc <- unlist(lapply(perf,function(x) x$auroc))
idx <- which.min(abs(auroc-mean(auroc)))
print(idx)



