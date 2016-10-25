#' collect results of brca simulations
#'
rm(list=ls())

resDir <- "~/tmp/TCGA_BRCA_runMany"
#dirList <- dir(path=resDir,pattern="R")
dirList <- paste("R",1:200,sep="")
dirList <- setdiff(dirList,c("R29","R134","R135"))

outmat <- matrix(nrow=length(dirList),ncol=6)
colnames(outmat) <- c("tp","fp","tn","fn","acc","ppv")
rownames(outmat) <- dirList

ctr <- 1
for (d in dirList){
        print(d)
        load(sprintf("%s/%s/FinalResults.Rdata",resDir,d))
        outmat[ctr,] <- unlist(out$testRes$perfStats)
        ctr <- ctr+1
}
outFile <- sprintf("%s/VM_results.txt",resDir)
write.table(outmat,file=outFile,sep="\t",col=T,row=T,quote=F)

